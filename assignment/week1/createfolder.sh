
# This script when run creates a template project folder structure at the current location.
# The folder contains seven directories: bash/ bin/ scripts/ doc/ data/ result/ assignment/


# Create project folder
mkdir ./project


# Create directories
for folders in bash bin scripts doc data
do
	mkdir ./project/$folders
done

for folders in result assignment
do
	mkdir ./project/$folders
	for week in week1 week2 week3 week4 week5 week6
	do
		mkdir ./project/$folders/$week
	done
done


# Create runall.sh as the the main driver script
# Create showfile.sh for showing all the directories and files in the project folder
# Store runall.sh showfile.sh in bash/
echo '# This script runs everything, from processing data, intermediate steps and final results.' >> ./project/bash/runall.sh
echo '# This script when run shows the directories and files in project folder.' >> ./project/bash/showfile.sh
echo '# Running this script needs to first locate to the project root directory' >> ./project/bash/showfile.sh
echo 'find .' >> ./project/bash/showfile.sh


# END
