import sys
#try to find module
if len(sys.argv) != 2:
    sys.exit("Please enter path to brex folder.")
sys.path.append(sys.argv[1])
print(sys.path)
try:
    import brex
except ImportError as e:
    sys.exit("Module brex not found. Make sure you entered the correct path.")
print("Loaded brex from Pythonpath")
try:
    import bamread
except ImportError as e:
    sys.exit("Module bamread not found. Make sure it is installed where python can find it.")
print("Loaded bamread from Pythonpath")
