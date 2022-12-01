cd
python -m venv mt_code
chmod +w ~/mt_code/bin/activate
# then we append the module load statement we have been using earlier to the activation script
#      NOTE:                                                      the -s will prevent the module load statement from cluttering your terminal or log file
echo module load lang/SciPy-bundle >> ~/mt_code/bin/activate

# in the end we need to protect the activate script from accidental modification:
chmod -w ~/mt_code/bin/activate

source ~/mt_code/bin/activate

pip install -r requirements.txt