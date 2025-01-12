# the script to create a pre-built file for the sample description

python3 scripts/make_sample_description.py -i input/peaks -o input/FILL_ME_UP.sample_description.csv

echo "
YOUR NEXT STEPS:
1) go to 'input' directory and find the created 'FILL_ME_UP.sample_description.csv' file
2) fill up empty fields in the file (replicas with the same 'Group name' will be processed together and saved in the individual output file)
3) rename it to 'sample_description.csv'
4) return to the main directory and run 'MetaPept.sh' script
"

