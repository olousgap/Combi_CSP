import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CombiCSP", # Replace with your own username
    version="0.1.0",
    author="N. Papadakis, G. Arnaoutakis",
    author_email="npapnet@gmail.com, g.e.arnaoutakis@gmail.com",
    description="A package for Concentrated Solar Collectors and Solar Tower Calculations ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)