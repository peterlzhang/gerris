#!/bin/sh 
# Generated automatically. Please modify crossref.sh.in.

path="/usr/local/share/gerris"

usage()
{
	cat <<EOF
Usage: crossref.sh [OPTIONS] FILE1 FILE2...

Creates cross-references

Options:
	[--url=URL] reference URL
        [--help]    displays this message and exits
EOF
	exit $1
}

if test $# -lt 1; then
	usage 1 1>&2
fi

while test $# -gt 1; do
  case "$1" in
  -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
  *) optarg= ;;
  esac

  case $1 in
    --url=*)
      url=$optarg
      ;;
    --help)
      usage 0 1>&2
      ;;
    --*)
      usage 0 1>&2
      ;;
      *)
      break
      ;;
  esac
  shift
done

keywords=`awk '{if ($1 == "gfs_keyword" && substr ($3,1,4) == "\"Gfs") print substr($3,2,length($3)-2); }' < $path/gfs.lang`

if test -d references; then :
else
    mkdir references
fi

for k in $keywords; do
    rm -f references/$k.html
    for f in $*; do
	gfsxref --url="$url/$f.html" $k < $f/$f.gfs >> references/$k.html
	cd $f
	for d in *; do
	    if test -d $d; then
		gfsxref --url="$url/$f.html#$d" $k < $d/$d.gfs >> ../references/$k.html
	    fi
	done
	cd ..
    done
done
