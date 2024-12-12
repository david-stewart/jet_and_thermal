#run runJetscape

echo "Number of command line parameters: $#"

# Check if there is at least one command line parameter
if [ $# -lt 1 ]; then
    echo "Error: Please provide a directory as a command line parameter."
    exit 1
fi

# Check if the directory exists
directory="$1"
if [ ! -d "$directory" ]; then
    echo "Error: The specified directory does not exist."
    exit 1
fi

# Continue with the rest of the script
ls -l $directory

for xml_file in `ls -1 $directory`; do
    in_file=$directory/$xml_file
    echo $in_file
    ls $in_file
    echo "hi $in_file"
    touch ${in_file}_CURRENTLY_RUNNING
    ./runJetscape $in_file
    rm ${in_file}_CURRENTLY_RUNNING
done
