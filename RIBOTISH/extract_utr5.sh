#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=4:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ribotish_prog
#SBATCH --output=ribotish_prog_%j.out
#SBATCH --error=ribotish_prog_%j.err

cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish
mkdir -p progressive_test

# Base command that we know works
BASE_CMD="ribotish quality \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
    -g /global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"

# Test 1: Base + output file
echo "Test 1: Base command with output file"
$BASE_CMD -o progressive_test/test1.txt
echo "Test 1 complete"

# Test 2: Add PDF output
echo -e "\nTest 2: Adding PDF output"
$BASE_CMD -o progressive_test/test2.txt -f progressive_test/test2.pdf
echo "Test 2 complete"

# Test 3: Add read length range
echo -e "\nTest 3: Adding read length range"
$BASE_CMD -o progressive_test/test3.txt -f progressive_test/test3.pdf -l 25,35
echo "Test 3 complete"

# Test 4: Add distance parameter
echo -e "\nTest 4: Adding distance parameter"
$BASE_CMD -o progressive_test/test4.txt -f progressive_test/test4.pdf -l 25,35 -d 40
echo "Test 4 complete"

# Test 5: Add parameter file output
echo -e "\nTest 5: Adding parameter file output"
$BASE_CMD -o progressive_test/test5.txt -f progressive_test/test5.pdf -l 25,35 -d 40 -r progressive_test/test5.para.py
echo "Test 5 complete"

# Test 6: Add number of processors
echo -e "\nTest 6: Adding processor count"
$BASE_CMD -o progressive_test/test6.txt -f progressive_test/test6.pdf -l 25,35 -d 40 -r progressive_test/test6.para.py -p 6
echo "Test 6 complete"

echo -e "\nAll tests complete. Check progressive_test directory for outputs."

# List all output files
echo -e "\nGenerated files:"
ls -l progressive_test/
