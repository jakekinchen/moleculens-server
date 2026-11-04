# OVH operations

Runtime: `/opt/moleculens`, project `moleculens-vps1`. Use the immutable image IDs in the current launch receipt. The private `.env` is root-owned, mode 0600. Only Caddy publishes ports. Keep legacy artifacts read-only and separate from new attempt storage.

`backup-database.sh` and the two systemd units are the deployed daily logical backup. Install the script at `/opt/moleculens/backup-database.sh` and the units under `/etc/systemd/system/`; then enable `moleculens-database-backup.timer`. Dumps retain seven versions under `/var/backups/moleculens`. Copy a matching dump and artifact archive off-host before a destructive maintenance operation. Provider automated backups are separate; a provider snapshot restore has not been demonstrated.

Run a logical restore check while job rows are stable (quiesce admission and wait for active attempts to finish first):

```sh
sudo python3 verify_restore.py --dump /var/backups/moleculens/<dump>.dump --output /opt/moleculens/receipts/restore-check.json --legacy-manifest /home/ubuntu/migration-backup/cache-initial-manifest.json
```

This creates and removes a uniquely named check database. It checks every Go job row, both historical tables, every referenced completed Go artifact, and optionally every historical artifact. It does not overwrite the live database or artifact volumes.

The current VPS profile has one chemistry worker and one conformer/naming worker. Each requires its configured attempt disk budget plus 2 GiB free before claiming work. API readiness/admission requires 2 GiB. Outer attempts are bounded at two hours, naming at eight seconds. Do not add replicas without revisiting total storage, CPU and RAM reserves. A syntactically accepted molecule is not a guarantee that its calculation fits these bounds. Check disk space, container memory/OOM events, queue state, and the latest backup/timer status when investigating a stalled queue.

`moleculens gc` defaults to a dry run. `gc --apply` removes only eligible unreferenced attempts older than 24 hours; completed results and uncertain publications remain retained. Do not change ages or remove referenced data to admit work. See the release runbook for rollback using the retained compatible Go API in quiesced mode.
