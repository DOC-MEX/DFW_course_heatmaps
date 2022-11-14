import requests
import json


server_url = "http://localhost:2000/grassroots/public_backend"
#server_url = "https://grassroots.tools/public_backend"

'''
Get study using id
returns JSON from backend
'''

def get_plot(id):
    plot_request = {
            "services": [{
                "so:name": "Search Field Trials",
                "start_service": True,
                "parameter_set": {
                    "level": "advanced",
                    "parameters": [{
                        "param": "ST Id",
                        "current_value": id
                    }, {
                        "param": "Get all Plots for Study",
                        "current_value": True
                    }, {
                        "param": "ST Search Studies",
                        "current_value": True
                    }]
                }
            }]
        }
    res = requests.post(server_url, data=json.dumps(plot_request))
    return json.dumps(res.json())

