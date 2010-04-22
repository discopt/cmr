#!/bin/sh
rm -rf autom4te.cache
libtoolize --force
aclocal -I m4
autoheader
automake --add-missing
autoconf
