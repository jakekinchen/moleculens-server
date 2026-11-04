package store

import (
	"context"
	"encoding/binary"
	"io"
	"net"
	"net/url"
	"strconv"
	"strings"
	"sync/atomic"
	"testing"
	"time"
)

// This proxy forwards the completion to a real PostgreSQL server, then drops
// the ReadyForQuery acknowledgement after the autocommit has completed. New
// connections pass normally, allowing the production readback path to resolve it.
func TestLostCommitAcknowledgementReadsBackRealPostgres(t *testing.T) {
	direct := testStore(t)
	ctx := context.Background()
	j, _, e := direct.Admit(ctx, "orbital", "ack-key", []byte(`{}`), 10)
	if e != nil {
		t.Fatal(e)
	}
	j, e = direct.Claim(ctx, "worker", time.Minute)
	if e != nil {
		t.Fatal(e)
	}
	listener, e := net.Listen("tcp", "127.0.0.1:0")
	if e != nil {
		t.Fatal(e)
	}
	defer listener.Close()
	config := direct.Pool.Config().ConnConfig
	backend := net.JoinHostPort(config.Host, strconv.Itoa(int(config.Port)))
	var dropped atomic.Bool
	go func() {
		for {
			front, e := listener.Accept()
			if e != nil {
				return
			}
			go func() {
				defer front.Close()
				back, e := net.Dial("tcp", backend)
				if e != nil {
					return
				}
				defer back.Close()
				var armed atomic.Bool
				go func() {
					defer back.Close()
					buf := make([]byte, 64<<10)
					tail := ""
					for {
						n, e := front.Read(buf)
						if n > 0 {
							tail += string(buf[:n])
							if strings.Contains(tail, "UPDATE moleculens_v2.jobs SET state='done'") && dropped.CompareAndSwap(false, true) {
								armed.Store(true)
							}
							if len(tail) > 8192 {
								tail = tail[len(tail)-8192:]
							}
							if _, err := back.Write(buf[:n]); err != nil {
								return
							}
						}
						if e != nil {
							return
						}
					}
				}()
				completionSeen := false
				for {
					header := make([]byte, 5)
					if _, e = io.ReadFull(back, header); e != nil {
						return
					}
					length := int(binary.BigEndian.Uint32(header[1:]))
					if length < 4 || length > 16<<20 {
						return
					}
					payload := make([]byte, length-4)
					if _, e = io.ReadFull(back, payload); e != nil {
						return
					}
					if armed.Load() && header[0] == 'C' && strings.HasPrefix(string(payload), "UPDATE 1") {
						completionSeen = true
					}
					if completionSeen {
						if header[0] == 'Z' {
							return
						}
						continue
					}
					if _, e = front.Write(append(header, payload...)); e != nil {
						return
					}
				}
			}()
		}
	}()
	u := url.URL{Scheme: "postgresql", Host: listener.Addr().String(), Path: "/" + config.Database, User: url.UserPassword(config.User, config.Password), RawQuery: "sslmode=disable"}
	proxied, e := Open(ctx, u.String())
	if e != nil {
		t.Fatal(e)
	}
	defer proxied.Pool.Close()
	deadline, cancel := context.WithTimeout(ctx, 10*time.Second)
	defer cancel()
	if e = proxied.Complete(deadline, j, []byte(`{"committed":true}`)); e != nil {
		t.Fatal(e)
	}
	if !dropped.Load() {
		t.Fatal("completion acknowledgement was not dropped")
	}
	actual, e := direct.Get(ctx, j.ID, j.Kind)
	if e != nil || actual.State != "done" {
		t.Fatal(actual, e)
	}
}
