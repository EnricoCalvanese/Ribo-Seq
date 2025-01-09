#!/usr/bin/env python3

# Save this as test_ribotish.py
import sys
import os
from importlib import import_module

# Try to import ribotish modules
try:
    ribo = import_module('ribotish.zbio.ribo')
    print("Successfully imported ribotish")
except ImportError as e:
    print(f"Could not import ribotish: {e}")
    sys.exit(1)

# Print the actual module location
print(f"Ribotish module location: {ribo.__file__}")

# Try to create a lenDis object with different distance formats
test_cases = [
    "-40,20",
    "[-40,20]",
    "-40 20",
    "-40:20"
]

for test in test_cases:
    print(f"\nTesting distance format: {test}")
    try:
        # Try to parse the distance string
        print(f"  Attempting to parse: {test}")
        # Print how ribotish sees this input
        print(f"  Raw input type: {type(test)}")
        print(f"  Raw input value: {test}")
    except Exception as e:
        print(f"  Error: {str(e)}")
