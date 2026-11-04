package engine

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"time"

	"github.com/jakekinchen/moleculens-server/internal/store"
)

var hashName = regexp.MustCompile(`^[0-9a-f]{64}$`)
var tokenName = regexp.MustCompile(`^[0-9a-f]{32}$`)

// CollectOrphans retains every completed result, every currently owned attempt,
// and ambiguous publications. It only removes retired, unreferenced attempts
// after 24h (well beyond the maximum 2h process deadline), up to 100 per run.
func CollectOrphans(ctx context.Context, s *store.Store, root string, apply bool) ([]string, error) {
	rows, e := s.Pool.Query(ctx, `SELECT cache_key,token FROM moleculens_v2.jobs WHERE token IS NOT NULL AND state IN ('pending','running','done')`)
	if e != nil {
		return nil, e
	}
	defer rows.Close()
	held := map[string]bool{}
	for rows.Next() {
		var key, token string
		if e = rows.Scan(&key, &token); e != nil {
			return nil, e
		}
		held[key+"/"+token] = true
	}
	if e = rows.Err(); e != nil {
		return nil, e
	}
	entries, e := os.ReadDir(root)
	if e != nil {
		return nil, e
	}
	candidates := []string{}
	for _, key := range entries {
		if !key.IsDir() || !hashName.MatchString(key.Name()) {
			continue
		}
		attempts, e := os.ReadDir(filepath.Join(root, key.Name()))
		if e != nil {
			return candidates, e
		}
		for _, attempt := range attempts {
			if !attempt.IsDir() || !tokenName.MatchString(attempt.Name()) {
				continue
			}
			rel := key.Name() + "/" + attempt.Name()
			if held[rel] {
				continue
			}
			path := filepath.Join(root, rel)
			info, e := attempt.Info()
			if e != nil {
				return candidates, e
			}
			if time.Since(info.ModTime()) < 24*time.Hour {
				continue
			}
			if _, e = os.Stat(filepath.Join(path, "publication-uncertain")); e == nil {
				continue
			} else if !os.IsNotExist(e) {
				return candidates, e
			}
			candidates = append(candidates, rel)
			if apply { // Recheck authority immediately before deletion; ambiguity aborts.
				var n int
				e = s.Pool.QueryRow(ctx, `SELECT count(*) FROM moleculens_v2.jobs WHERE cache_key=$1 AND token=$2 AND state IN ('pending','running','done')`, key.Name(), attempt.Name()).Scan(&n)
				if e != nil {
					return candidates, e
				}
				if n > 0 {
					return candidates, fmt.Errorf("attempt acquired a reference; aborting collection")
				}
				if e = os.RemoveAll(path); e != nil {
					return candidates, e
				}
			}
			if len(candidates) >= 100 {
				return candidates, nil
			}
		}
	}
	return candidates, nil
}
