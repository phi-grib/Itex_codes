# Itex_codes

Python scripts that generate and update the CII database from Inditex project.

Installing itex_codes:<br>
* $ conda env create -f environment.yml<br>
* $ conda activate inditex<br>
* $ python setup.py install<br>
* $ pip install -r requirements.txt<br>

Once this is done, you should be able to run the code without any warning in the exports.

*Note:* the code under CreateDB it's not used anymore, since it was written for a one time use and after CII v0 was released only three functions under functions_to_process.py are still used in Update_CII.py
