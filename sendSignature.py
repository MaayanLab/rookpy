

import requests
import json


payload = {}
payload["direction"] = "similar"
payload["signatureName"] = "Example_similar"
payload["species"] = "human"
payload["type"] = "geneset"
payload["upgenes"] =  ["AARS", "ACLY", "KIAA0907", "KDM5A", "CDC25A", "EGR1", "GADD45B", "RELB", "TERF2IP", "SMNDC1"]
payload["downgenes"] =  ["ACTN1", "ACTG1", "SCCPDH", "KIF20A", "FZD7", "USP22", "PIP4K2B", "CRYZ", "GNB5", "EIF4EBP1", "PHGDH"]


upload_url = 'https://maayanlab.edu/rookpy/signature'

upload_url = "http://localhost:5055/rookpy/signatures"

print(json.dumps(payload))

r = requests.post(upload_url, json.dumps(payload))

print(r.text)




































