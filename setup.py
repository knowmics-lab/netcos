# drug_repurposing/setup.py
from setuptools import setup, find_packages

def read_requirements():
    with open('requirements.txt', 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name="drug_repurposing",
    version="0.1.0",
    author="L-F-S",
    description="Network-Enhanced Connectivity Score for drug repurposing",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.7",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "jupyter",
            "ipykernel",
            "pytest",
        ]
    }
)