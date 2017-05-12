#!/bin/sh

# BEAST executable
BEAST=beast

cd prior_sampling_nonstarbeast
for xml in nonstarbeast_3taxon.xml nonstarbeast_4taxon.xml; do
    $BEAST -overwrite $xml
done
cd ..

cd prior_sampling_starbeast
for xml in starbeast_3taxon.xml starbeast_4taxon.xml; do
    $BEAST -overwrite $xml
done
cd ..
