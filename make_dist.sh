#!/bin/bash

VERSION=1.2b
NAME=unimodularity-library-Polymake-${VERSION}

rm -rf ${NAME}
mkdir ${NAME}
cp -r URI LICENSE_1_0.txt AUTHORS ChangeLog NEWS README apps/ examples/ ${NAME}/
tar czf ${NAME}.tar.gz ${NAME}
rm -rf ${NAME}
