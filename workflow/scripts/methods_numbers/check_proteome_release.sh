#!/bin/bash
# Which Ensembl release are the SAMap proteomes from? No download record exists, so compare the local files'
# modification times (preserved from the FTP server) with the Ensembl FTP listings of candidate releases.
# Read-only; prints to stdout. Usage: bash check_proteome_release.sh [releases...]   (default: 111 112 113 114)
PF=/n/lab_storage/eddy_lab/users/nhilgert/data/proteomes/fasta
echo "## local (mtime in UTC)"
for f in "$PF"/Homo_sapiens.GRCh38.pep.all.fa "$PF"/Mus_musculus.GRCm39.pep.all.fa; do
  echo "$(TZ=UTC date -r "$f" '+%Y-%m-%d %H:%M') $(basename "$f")"
done
for rel in ${@:-111 112 113 114}; do
  for sp in homo_sapiens mus_musculus; do
    # listing times are UK local time (Europe/London)
    line=$(timeout 12 curl -s "https://ftp.ensembl.org/pub/release-$rel/fasta/$sp/pep/" | grep -E 'pep\.all\.fa\.gz' | sed 's/<[^>]*>/ /g' | tr -s ' ' | head -1)
    set -- $line
    if [ -n "$2" ]; then
      utc=$(TZ=UTC date -d "TZ=\"Europe/London\" $2 $3" '+%Y-%m-%d %H:%M')
      echo "## release $rel $sp: $utc UTC ($1, $4)"
    else
      echo "## release $rel $sp: listing unavailable"
    fi
  done
done
