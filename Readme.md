# for-lab-needs

This is a repository to address minor bioinformatics-related needs of Systems Biology lab at IMBG of NASU.

## Features 

	* some.py - splits double strain DNA into single strain semi-overlapping pieces of reasonably similar melting temperature in segments of overlap; this allows to recreate original DNA from smaller pieces by mixing them together and iteratively heating/cooling the mixture.


## Getting Started

	```
	python some.py ATATAGATTACAATATA
	```
	or just
	```
	python some.py
	```
	If unspecified as a parameter sequence shall be retreived from `input.txt` file of the same folder as `some.py`

### Prerequisites

python3, pip, virtualenv

### Installing

	```
	git clone https://github.com/Sashkow/for-lab-needs.git
	cd for-lab-needs
	virtualenv -p python3 env
	source env/bin/activate
	pip install -r requirements.txt
	python some.py
	```

<!-- ## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system
 --> 
## Built With

* [Python 3.5.2](https://www.python.org/downloads/release/python-352/) - Programming language


<!-- ## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 
 
## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.
-->
## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details


## Acknowledgments

	* Wiss for the idea and contribution
