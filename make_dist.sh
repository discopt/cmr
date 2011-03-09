#!/bin/bash

VERSION=0.9
NAME=unimodularity-library-polymake-${VERSION}

rm -rf ${NAME}
mkdir ${NAME}
cp -r URI LICENSE_1_0.txt AUTHORS ChangeLog NEWS README apps/ examples/ ${NAME}/
tar czf ${NAME}.tar.gz ${NAME}
rm -rf ${NAME}
