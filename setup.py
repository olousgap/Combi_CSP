import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    'matplotlib',
    'scipy',
    'ipykernel',
    'pandas',
    'pvlib',
    'iapws',
    'numpy-financial',
    'seaborn',
]

test_requirements = [
    'pytest',
    # 'pytest-pep8',
    # 'pytest-cov',
]


setuptools.setup(
    name="CombiCSP", # Replace with your own username
    version="1.5.0",
    author="N. Papadakis, G. Arnaoutakis",
    author_email="npapnet@gmail.com, g.e.arnaoutakis@gmail.com",
    description="A package for Concentrated Solar Collectors and Solar Tower Calculations ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/npapnet/Combi_CSP/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements,
    tests_require=test_requirements,
    python_requires='>=3.8',
)