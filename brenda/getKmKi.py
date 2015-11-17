#!/usr/bin/python
import string
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file
import re

#1) Usage with WSDL (for extracting the URL of the endpoint)
#wsdl = "http://www.brenda-enzymes.info/soap2/brenda.wsdl"
#client = WSDL.Proxy(wsdl)
## Km ##
# resultString = client.getKmValue("ecNumber*1.1.1.1#organism*Homo sapiens")
## Ki ##
# resultString = client.getKiValue("ecNumber*1.1.1.1#organism*Homo sapiens")

#2) Usage without WSDL
endpointURL = "http://www.brenda-enzymes.info/soap2/brenda_server.php"
client = SOAPProxy(endpointURL)

########
## Choose one of the lines in this block
## Km ##
#resultString = client.getKmValue("organism*Mus musculus")
resultString = client.getKmValue("organism*Homo sapiens")
#resultString = client.getKmValue("organism*Saccharomyces cerevisiae")
#resultString = client.getKmValue("organism*Escherichia coli")
#resultString = client.getKmValue("organism*Escherichia coli K-12")
#resultString = client.getKmValue("ecNumber*1.1.1.1#organism*Homo sapiens") #for specific rxn
## Ki ##
#resultString = client.getKiValue("organism*Mus musculus")
#resultString = client.getKiValue("organism*Homo sapiens")
#resultString = client.getKiValue("organism*Saccharomyces cerevisiae")
#resultString = client.getKiValue("organism*Escherichia coli")
#resultString = client.getKiValue("organism*Escherichia coli K-12")
#resultString = client.getKiValue("ecNumber*1.1.1.1#organism*Homo sapiens") #for specific rxn
########

# For further parsing
resultString=re.sub('!','\n',resultString)
resultString=re.sub('\*','#',resultString)

print resultString.encode('ascii','xmlcharrefreplace')