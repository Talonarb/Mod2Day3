import os
import sys
os.chdir(os.path.dirname(sys.argv[0]))
PROJECT_PATH = os.path.dirname(os.getcwd())
SOURCE_PATH = os.path.join(
    PROJECT_PATH,"src"
)
sys.path.append(SOURCE_PATH)