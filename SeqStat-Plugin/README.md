# Plugin-Visual-Test
A very small test plugin to test the visualisation interface is working for SynBioHub. Could be the basis for other python based visualisation plugins.

# Install
## Using docker
Run `docker run --publish 8080:5000 --detach --name python-test-plug synbiohub/plugin-visual-test:snapshot`
Check it is up using localhost:8080.  

## Using Python
Run `pip install -r requirements.txt` to install the requirements. Then run `FLASK_APP=app python -m flask run`. A flask module will run at localhost:5000/.

# Using SeqStat with a terminal
Bellow are some terminal commands wich can be used to access a webpage instance of a SynBioHub part and run it through the plugin. Note that my file location is 
/Users/benjacobsen/Desktop/Plugin-Visual-Test/output.html
That must be changed for any other instance of this plugin being ran. Also the plugin folder is now called "SeqStat-Plugin"

FLASK_APP=app

flask run

curl -X POST http://127.0.0.1:5000/run \
-H "Content-Type: application/json" \
-d '{
    "top_level": "http://example.org/top_level_component",
    "complete_sbol": "http://example.org/sbol/document.sbol",
    "instanceUrl": "http://example.org/sequence_instance",
    "size": 500,
    "type": "Sequence",
    "shallow_sbol": ""
}' > visualization.html

curl -X POST http://127.0.0.1:5000/run \
-H "Content-Type: application/json" \
-d '{
  "top_level": "https://synbiohub.org/public/igem/BBa_E0040/1",
  "complete_sbol": "https://synbiohub.org/public/igem/BBa_E0040/1/sbol",
  "instanceUrl": "https://synbiohub.org/",
  "size": 500,
  "type": "Component",
  "shallow_sbol": ""
}' > /Users/benjacobsen/Desktop/Plugin-Visual-Test/output.html

curl -X POST http://127.0.0.1:5000/run \
-H "Content-Type: application/json" \
-d '{
  "top_level": "https://synbiohub.org/public/Eco1C1G1T1/BBa_B0064_rbs/1",
  "complete_sbol": "https://synbiohub.org/public/Eco1C1G1T1/BBa_B0064_rbs/1/sbol",
  "instanceUrl": "https://synbiohub.org/",
  "size": 500,
  "type": "Component",
  "shallow_sbol": ""
}' > /Users/benjacobsen/Desktop/Plugin-Visual-Test/output.html

curl -X POST http://127.0.0.1:5000/run \
-H "Content-Type: application/json" \
-d '{
  "top_level": "https://synbiohub.org/public/igem/BBa_K1444011/1",
  "complete_sbol": "https://synbiohub.org/public/igem/BBa_K1444011/1/sbol",
  "instanceUrl": "https://synbiohub.org/",
  "size": 500,
  "type": "Component",
  "shallow_sbol": ""
}' > /Users/benjacobsen/Desktop/Plugin-Visual-Test/output.html

curl -X POST http://127.0.0.1:5000/run \
-H "Content-Type: application/json" \
-d '{
  "top_level": "https://synbiohub.org/public/igem/BBa_K1481006/1",
  "complete_sbol": "https://synbiohub.org/public/igem/BBa_K1481006/1/sbol",
  "instanceUrl": "https://synbiohub.org/",
  "size": 500,
  "type": "Component",
  "shallow_sbol": ""
}' > /Users/benjacobsen/Desktop/Plugin-Visual-Test/output.html
