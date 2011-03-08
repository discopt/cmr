#!/bin/sh

BOOST_URL=http://downloads.sourceforge.net/project/boost/boost/1.43.0/boost_1_43_0.tar.bz2
BOOST_TGZ=boost_1_43_0.tar.bz2

BOOST_FOLDERS=
for FOLDER in numeric/ublas serialization config mpl type_traits iterator tuple integer multi_index property_map concept smart_ptr optional archive dynamic_bitset parameter range bind python graph detail preprocessor utility exception unordered pending functional random logic; do
  BOOST_FOLDERS="${BOOST_FOLDERS} boost_1_43_0/boost/${FOLDER}"
done

cd src
rm -rf boost
wget ${BOOST_URL} -O ${BOOST_TGZ}
tar xjvf ${BOOST_TGZ} --wildcards --no-wildcards-match-slash 'boost_1_43_0/boost/*.hpp' ${BOOST_FOLDERS}
mv boost_1_43_0/boost boost
rm -rf boost_1_43_0
rm ${BOOST_TGZ}
