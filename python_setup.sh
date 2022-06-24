cd
python -m venv venv
chmod +w ~/venv/bin/activate
# then we append the module load statement we have been using earlier to the activation script
#      NOTE:                                                      the -s will prevent the module load statement from cluttering your terminal or log file
echo module load lang/SciPy-bundle >> ~/venv/bin/activate

# in the end we need to protect the activate script from accidental modification:
chmod -w ~/venv/bin/activate

pip install Bio
pip install pyranges
pip install gffpandas