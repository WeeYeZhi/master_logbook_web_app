#!/bin/bash

set -eux -o pipefail

echo "Launching LongStitch test runs..."
echo "Ensure all dependencies are in your PATH"

echo "Running LongStitch steps to scaffold the CPB genome assembly..."
longstitch -B run draft=test_scaffolds1 reads=test_reads1 G=560450000 longmap=ont t=48 out_prefix=CPB

set +x

echo ""
echo "Done tests! Compare your generated files with the files in the expected_outputs folder to ensure the tests were successful."
echo "Note: Files were generated using Tigmint v1.2.5, ntLink v1.3.3 and ARCS v1.2.4"
echo ""

echo "Final scaffold files found in: test_scaffolds1.k24.w100.tigmint-ntLink.longstitch-scaffolds.fa test_scaffolds2.k32.w150.tigmint-ntLink-arks.longstitch-scaffolds.fa"
echo "Tip: compare the fasta files using abyss-fac. Some files might appear slightly different using command line tools like 'cmp' due to sequences being output on different strands."

