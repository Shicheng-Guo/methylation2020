import urllib2
from urllib2 import Request, urlopen, URLError
import json
import math
import sys
import os

RunID = raw_input("Enter run ID:")
AccessToken = raw_input("Enter access token:")

def restrequest(rawrequest):
	#request = Request(rawrequest)
	#print rawrequest
	try:
		opener = urllib2.build_opener()
		opener.addheaders = [('User-agent', 'Mozilla/5.0')]
		response = opener.open(rawrequest)
		#response = urlopen(request)
        	json_string = response.read()
        	json_obj = json.loads(json_string)

	except URLError, e:
    		print 'Got an error code:', e
		sys.exit()

	return json_obj

def downloadrestrequest(rawrequest,path):
	dirname = os.path.dirname(path)
	basename = os.path.basename(path)
	if dirname != '':
		if not os.path.isdir(dirname):
			os.system('mkdir -p %s' %(dirname))
	
	request = (rawrequest)
	#print request
        try:
                opener = urllib2.build_opener()
                opener.addheaders = [('User-agent', 'Mozilla/5.0')]
                response = opener.open(rawrequest)

		outfile = open(path,'wb')
		outfile.write(response.read())
		outfile.close()
		

        except URLError, e:
                print 'Got an error code:', e
                sys.exit()


#AccessToken = '36c60e666fa24ca5bcf6ed39a5fd6d6a'
		
request = 'http://api.basespace.illumina.com/v1pre3/runs/%s/files?access_token=%s' %(RunID,AccessToken)

json_obj = restrequest(request)
totalCount = json_obj['Response']['TotalCount']

noffsets = int(math.ceil(float(totalCount)/1000.0))

hreflist = []
pathlist = []
filenamelist = [] 
for index in range(noffsets):
	offset = 1000*index
	request = 'http://api.basespace.illumina.com/v1pre3/runs/%s/files?access_token=%s&limit=1000&Offset=%s' %(RunID,AccessToken,offset) 
	json_obj = restrequest(request)
	nfiles = len(json_obj['Response']['Items'])
	for fileindex in range(nfiles):
		href = json_obj['Response']['Items'][fileindex]['Href']
		hreflist.append(href)
		path = json_obj['Response']['Items'][fileindex]['Path']
		pathlist.append(path)

for index in range(len(hreflist)):
	request = 'http://api.basespace.illumina.com/%s/content?access_token=%s'%(hreflist[index],AccessToken)
	print 'downloading %s' %(pathlist[index]) 
	downloadrestrequest(request, pathlist[index])


