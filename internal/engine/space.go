package engine

import (
	"fmt"
	"syscall"
)

// StorageReserve keeps the API and small jobs usable while a large worker owns
// its scratch allocation. Each worker separately requires its full attempt budget
// plus this reserve before claiming. The VPS profile runs one large worker and
// one 1-GiB worker; provision additional reserve for additional large workers.
const StorageReserve uint64 = 2 << 30

func CheckSpace(root string) error {
	return checkSpace(root, StorageReserve)
}

func checkSpace(root string, required uint64) error {
	var s syscall.Statfs_t
	if e := syscall.Statfs(root, &s); e != nil {
		return e
	}
	free := uint64(s.Bavail) * uint64(s.Bsize)
	if free < required {
		return fmt.Errorf("artifact filesystem has less than %d MiB free", required>>20)
	}
	return nil
}
