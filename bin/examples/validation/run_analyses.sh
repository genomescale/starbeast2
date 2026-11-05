#!/bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BEAST="$PROJECT_ROOT/bin/beast"

# bin/beast always cds to $PROJECT_ROOT before running (so Maven can resolve
# the module path), so BEAST's $(filebase) outputs land in $PROJECT_ROOT
# rather than next to the XML. Run with an absolute XML path, then move the
# generated files (named after the XML's basename) into the XML's directory.
run_analysis() {
    dir=$1
    xml=$2
    base=$(basename "$xml" .xml)

    "$BEAST" -overwrite "$SCRIPT_DIR/$dir/$xml"

    for f in "$PROJECT_ROOT/$base".*; do
        [ -e "$f" ] && mv "$f" "$SCRIPT_DIR/$dir/"
    done
}

for xml in nonstarbeast_3taxon.xml nonstarbeast_4taxon.xml; do
    run_analysis prior_sampling_nonstarbeast $xml
done

for xml in starbeast_3taxon.xml starbeast_4taxon.xml; do
    run_analysis prior_sampling_starbeast $xml
done
