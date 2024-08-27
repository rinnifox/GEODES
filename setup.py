import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geodes",
    version="0.0.1",
    author="Karina Pats",
    author_email="karina.m.pats@gmail.com",
    description="GEODES package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rinnifox/GEODES",
    project_urls={
        "Bug Tracker": "https://github.com/rinnifox/GEODES/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Ubuntu",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9",
)
