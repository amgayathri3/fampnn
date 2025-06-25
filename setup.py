from setuptools import setup, find_packages

setup(
    name='fampnn',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'fampnn-design = run_fampnn_wrapper:main'
        ],
    },
)
