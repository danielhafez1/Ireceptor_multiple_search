import urllib.request, urllib.parse
import requests
import pandas as pd
import numpy as np
import argparse
import json
import os, ssl
import sys
import time
import yaml
from argparse import ArgumentParser
from tqdm import tqdm
import pickle
import datetime
import urllib3
urllib3.disable_warnings()
requests.urllib3.disable_warnings()
print("\n Started at --> {0} \n".format(datetime.datetime.now()))

### Never see in i-receptor but i found one sequence with this URL while looking for "CAVSGSGGYQKVTF"
#"https://stats-staging.ireceptor.org"

def get_list_repositories():
    list_rep = ["https://roche-airr.ireceptor.org","https://agschwab.uni-muenster.de","http://covid19-1.ireceptor.org", "http://covid19-2.ireceptor.org", "http://covid19-3.ireceptor.org", "http://covid19-4.ireceptor.org","https://ipa1.ireceptor.org", "https://ipa2.ireceptor.org", "https://ipa3.ireceptor.org", "https://ipa4.ireceptor.org","http://ipa5-staging.ireceptor.org", "https://vdjserver.org", "https://scireptor.dkfz.de"]
    
    df = pd.DataFrame({"URL": list_rep})
    return df


def processQuery(query_url, header_dict, query_json="{}", verbose=False):
    # Build the required JSON data for the post request. The user
    # of the function provides both the header and the query data
    # Encode the JSON for the HTTP requqest
    query_json_encoded = query_json.encode('utf-8')
    # Try to connect the URL and get a response. On error return an
    # empty JSON array.
    try:
        # Build the request
        request = urllib.request.Request(query_url, query_json_encoded, header_dict)
        # Make the request and get a handle for the response.
        response = urllib.request.urlopen(request)
        # Read the response
        url_response = response.read()
        # If we have a charset for the response, decode using it, otherwise assume utf-8
        if not response.headers.get_content_charset() is None:
            url_response = url_response.decode(response.headers.get_content_charset())#print("########URL_RESPONSE",url_response)
        else:
            url_response = url_response.decode("utf-8")
    except urllib.error.HTTPError as e:
        print('ERROR: Server could not fullfil the request to ' + query_url)
        print('ERROR: Error code = ' + str(e.code) + ', Message = ', e.read())
        return json.loads('[]')
    except urllib.error.URLError as e:
        print('ERROR: Failed to reach the server')
        print('ERROR: Reason =', e.reason)
        return json.loads('[]')
    except Exception as e:
        print('ERROR: Unable to process response')
        print('ERROR: Reason =' + str(e))
        return json.loads('[]')
    try:
        json_data = json.loads(url_response)
    except json.decoder.JSONDecodeError as error:
        if force:
            print("WARNING: Unable to process JSON response: " + str(error))
            if verbose:
                print("Warning: URL response = " + url_response)
            return json.loads('[]')
        else:
            print("ERROR: Unable to process JSON response: " + str(error))
            if verbose:
                print("ERROR: URL response = " + url_response)
            return json.loads('[]')
    except Exception as error:
        print("ERROR: Unable to process JSON response: " + str(error))
        if verbose:
            print("ERROR: JSON = " + url_response)
        return json.loads('[]')
    # Return the JSON data
    return json_data

def getHeaderDict():
    # Set up the header for the post request.
    header_dict = {'accept': 'application/json',
                   'Content-Type': 'application/json'}
    return header_dict

def initHTTP():
    # Deafult OS do not have create cient certificate bundles. It is
    # easiest for us to ignore HTTPS certificate errors in this case.
    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
        getattr(ssl, '_create_unverified_context', None)):
        ssl._create_default_https_context = ssl._create_unverified_context

def generateQuery(field, value):
    query_str = "{ \"filters\": { \"op\":\"contains\", \"content\": { \"field\":\"%s\", \"value\":\"%s\" } } }"%(field, value)
    #query_str = "{ \"filters\": { \"op\":\"contains\", \"content\": { \"field\":\"locus",\"value\":\"TRB\"} },{\"op\":\"contains\", \"content\": { \"field\":\"%s",\"value\":\"%s\"}}}"%(field,value)
    return query_str


def generateQuery_VDJ(field, value,first):
    query_str = "{ \"filters\": { \"op\":\"contains\", \"content\": { \"field\":\"%s\", \"value\":\"%s\" } },\"from\": %s,\"size\": 1000 }"%(field, value,first)
    #query_str = "{ \"filters\": { \"op\":\"contains\", \"content\": { \"field\":\"locus",\"value\":\"TRB\"} },{\"op\":\"contains\", \"content\": { \"field\":\"%s",\"value\":\"%s\"}}}"%(field,value)
    return query_str

def searchCDR3_single_seq(url, cdr3_seq, verbose):
    # Ensure our HTTP set up has been done.
    initHTTP()
    # Get the HTTP header information (in the form of a dictionary)
    header_dict = getHeaderDict()
    # Open the file that contains the list of CDR3s to search
    # Get the CDR3 list from the column header.
    # Build the full URL combining the URL and the entry point.
    query_url = url
    # Iterate over the CDR3s
    #df_info_total =pd.DataFrame()
    df_data_total =pd.DataFrame()
    #print("AAAAAAABBBBBBBB",cdr3_row[cdr3_header]);print(index)
    cdr3_query = generateQuery("junction_aa", cdr3_seq)
    if verbose:
        print('INFO: Performing query: ' + str(cdr3_query))
    # Perform the query.
    try:
        query_json = processQuery(query_url, header_dict, cdr3_query, verbose)
    except:
        print("\n####### NO seq for {0}. {1}".format(cdr3_seq, query_json))
        return 1

    json_data = query_json['Rearrangement']
    
    df_data = pd.json_normalize(json_data)
    df_data_total = df_data_total.append(df_data)
    
    if query_url == "https://vdjserver.org/airr/v1/rearrangement" and len(df_data_total) == 1000:
    	z = 1000
    	while len(json_data) == 1000:
    		cdr3_query_VDJ = generateQuery_VDJ("junction_aa", cdr3_seq,first = z)
    		
    		try:
    			query_json_VDJ = processQuery(query_url, header_dict, cdr3_query_VDJ, verbose)


    		except:
    			print("ERROR ONE THE WHILE LOOP FOR VDJ SERVER")
    			break
    		json_data = query_json_VDJ['Rearrangement']
    		
    		df_data_VDJ = pd.json_normalize(json_data)
    		
    		df_data_total = df_data_total.append(df_data_VDJ)
    		
    		z = z + 1000

    
    #if "single_cell" in ID_study:
    #data["Single_cell"] = ID_study["single_cell"]
    #df_info_total = df_info_total.append(df_info)
    if not df_data_total.empty :
        ##if "cell_id" in df_data_total.columns:
       	    	##data = pd.DataFrame(df_data_total[["junction_aa","repertoire_id","data_processing_id","locus","cell_id","sequence","sequence_id"]])###"D GENE DOESNT WORK FOR VDJ SERVER"
 

        ##if 'v_gene' in df_data_total.columns and 'j_gene' in df_data_total.columns:
        	#data = pd.DataFrame(df_data_total[["junction_aa","repertoire_id","data_processing_id","locus","sequence","sequence_id","v_gene","j_gene"]])###"D GENE DOESNT WORK FOR VDJ SERVER"
        
        ##if 'v_gene' in df_data_total.columns and 'j_gene' in df_data_total.columns and 'cell_id' in df_data_total.columns:

        	##data = pd.DataFrame(df_data_total[["junction_aa","repertoire_id","data_processing_id","locus","cell_id","sequence","sequence_id","v_gene","j_gene"]])###"D GENE DOESNT WORK FOR VDJ SERVER"

        data = df_data_total
        data["URL"] = query_url
        with pd.option_context('display.max_rows', None,'display.max_columns', None,'display.precision', 3,):
            return data
    

    else:
        data = pd.DataFrame()
        return data


def getRepertoires(repertoire_url):
    initHTTP()
    url = repertoire_url
    url_response = requests.post(url,verify = False)

    #repertoire_info = dict()
    json_data= url_response.json()
    repertoire_info = json_data["Repertoire"]
    #print("Searching in {0} repertoires ".format(len(repertoire_info)))
    ID_study = pd.json_normalize(repertoire_info)
    #data = pd.DataFrame(ID_study[["repertoire_id","subject.subject_id","study.study_id","subject.diagnosis","sample"]])
    data = pd.DataFrame(ID_study)

    #MHC COLUMNS: N.B different name for mhc columns in different repositories.... -> adapt them
    mhc_col = "MHC"
    if 'vdjserver' in repertoire_url:
        mhc_col = 'subject.mhc'
    if 'airr-covid-19' in repertoire_url:
        mhc_col = 'subject.genotype.mhc_genotype_set.mhc_genotype_set_id'
    if mhc_col in ID_study.columns:
        data[mhc_col] = ID_study[mhc_col]

    #print("##########################LENGTH",len(repertoire_info))
    #Diagnosis = []
    #for repertoire in repertoire_info:
        #ID_study = ID_study.append(I)
        #for i in repertoire["subject"]["diagnosis"]:#["disease_diagnosis"]:
            #Diagnosis.append(i["disease_diagnosis"]["label"])
            #print(i["disease_diagnosis"]["label"])
            #for j in i:
                #print(j)
            #for j in i["disease_diagnosis"]:
                #print(j.get('label'))
                #Diagnosis.append(i["disease_diagnosis"]["label"])
    #data["DIAGNOSIS"] = Diagnosis
    with pd.option_context('display.max_rows', None,'display.max_columns', None,'display.precision', 3,):
        return data
    #repertoire_ids = []
    #keys = ['study']


def getArguments():
    # Set up the command line parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=""
    )
    # The URL for the repository to test
    #parser.add_argument("repository_url_file")
    # The API entry point to use
    parser.add_argument("cdr3_file")
    # Comma separated list of query files to test.
    # Verbosity flag
    parser.add_argument(
        "--repertoire_api",
        dest="repertoire_api",
        default="/airr/v1/repertoire",
        help="The repertoire API entry point. Defaults to '/airr/v1/repertoire'/")
    parser.add_argument(
        "--rearrangement_api",
        dest="rearrangement_api",
        default="/airr/v1/rearrangement",
        help="The Rearrangement API entry point. Defaults to '/airr/v1/rearrangement'/")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Run the program in verbose mode.")
    # Parse the command line arguements.
    options = parser.parse_args()
    return options

if __name__ == "__main__":

    # Parser
    parser = ArgumentParser()
    #file with CDR3 data
    parser.add_argument("--cdr3_file",  default = './test.csv')
    #only cdr3 with the same V,J gene than the template TCR
    parser.add_argument( "--same_VJ", action="store_true", help="keep only seq with the same V,J genes (perfect match)")
    #path out folder
    parser.add_argument("--out_folder",  default = 'results')
    #model params
    parser.add_argument( "--repertoire_api", dest="repertoire_api", default="/airr/v1/repertoire", help="The repertoire API entry point. Defaults to '/airr/v1/repertoire'/")
    parser.add_argument( "--rearrangement_api", dest="rearrangement_api", default="/airr/v1/rearrangement", help="The Rearrangement API entry point. Defaults to '/airr/v1/rearrangement'/")
    parser.add_argument( "-v", "--verbose", action="store_true", help="Run the program in verbose mode.")
    options = parser.parse_args()

    if os.path.isdir(options.out_folder):
        print("ERROR! Output folder {0} already exists! Please choose another name".format(options.out_folder))
        sys.exit(1)

    #get the list of all repositories url
    df_all_repositories = get_list_repositories()

    dfs = []
    repertoire_entry_point = "/airr/v1/repertoire"
    rearrangement_entry_point = "/airr/v1/rearrangement"

    print("### Loading cdr3_file {0}. Starting search in iReceptor ####".format(options.cdr3_file))
    #print("N.B check format V,J genes in your dataset. iReceptor uses IMGT nomenclature \n(e.g. TRAV-8*04 OK || TRAV08*04 or TRAV8-2/4-8 not OK)")
    
    os.makedirs(options.out_folder)
    test_name_out = os.path.join(options.out_folder, options.cdr3_file.split("./")[-1].replace(".csv", "_results.pickle"))




    df_cdr3 = pd.read_csv(options.cdr3_file)
    
    for chain in ['TRA', 'TRB']:
        df_chain = df_cdr3.dropna(subset = ['cdr3_{0}'.format(chain)])
        for _, df_cdr3_tmp in df_chain.iterrows():
            d_chain_seq = {}
            df_merged = pd.DataFrame()
            df_final = pd.DataFrame()
            cdr3 = df_cdr3_tmp['cdr3_{0}'.format(chain)]
            v = df_cdr3_tmp['{0}V'.format(chain)]
            j = df_cdr3_tmp['{0}J'.format(chain)]
            print("# chain={0}\t{1}\t{2}\t{3}".format(chain, v, cdr3, j))
            if(len(cdr3) < 5):
                print(" Not valid seq, skipping it...")
                continue
            all_row_rep = df_all_repositories['URL'].values
            #print(all_row_rep)
            pbar = tqdm(all_row_rep)
            for row_rep in pbar:
                try:
                    pbar.set_description("Processing {0}".format(row_rep))
                    df_repertoires = getRepertoires(row_rep + repertoire_entry_point) #options.repertoire_api)
                except:
                    print('{0} repository will be passed try again later (unreachable)'.format(row_rep))
                    continue
                df_cdr3_search  = searchCDR3_single_seq(row_rep + rearrangement_entry_point, cdr3, False) #options.verbose)
                if len(df_cdr3_search) == 0:
                    #print("%%% NO MATCHES in this URL: {0}".format(row_rep))
                    continue
                df_merged  = pd.merge(df_repertoires, df_cdr3_search ,on="repertoire_id")
                ###NEED TO DISCUSS IF I REMOVE DUPLETS
                #df_merged = df_merged.drop_duplicates(subset = "sequence_id")#,keep = False, inplace = True)
                #print("###### FOUND {0} sequences in {1}".format(df_merged.shape[0], row_rep))
                if options.same_VJ:
                    v_gene = df_cdr3_tmp['{0}V'.format(chain)]
                    j_gene = df_cdr3_tmp['{0}J'.format(chain)]
                    #print("Considering only TCR with the same V ({0}) and J ({1}) genes".format(v_gene, j_gene))
                    df_merged = df_merged.loc[df_merged['v_gene'].str.contains(v_gene)]
                    if len(df_merged) > 0:
                        df_merged = df_merged.loc[df_merged['j_gene'].str.contains(j_gene)]
                    #print("Keeping {0} sequences".format(len(df_merged)))
                if len(df_merged) == 0:
                    #print("FOUND NO MATCHES")
                    continue
                if row_rep=='https://vdjserver.org' and len(df_merged) == 1000:
                    print("#### BE CERFULL AS THERE COULD BE MISSING SEQUENCES FROM THE VDJ SERVER FOR {0}".format(cdr3))
                df_final = df_final.append(df_merged)

                #except:
                #    print("ERROR with {0}".format(row_rep))
                #    continue
            
            df_cdr3_irec = df_final
            df_cdr3_irec = df_cdr3_irec.reset_index(drop = True)
            print("\t{0} sequences were found".format(len(df_cdr3_irec)))
            #add to final dictionary
            key = (chain, cdr3)
            if not df_cdr3_irec.empty:
                d_chain_seq[key] = df_cdr3_irec


            with open(test_name_out, 'wb') as filehandler:
            	pickle.dump(d_chain_seq, filehandler)

            with open(test_name_out, 'rb') as f:
            	model = pickle.load(f)

            print("+++ pickle object saved and reloaded, saving the .csv file for each TCRs, everything seems OK! +++++")

            for key, df_irec in model.items():
            	name_df = 'data_' + '_'.join(key) + '.csv'
            	path_df = os.path.join(options.out_folder, name_df)
            	model[key].to_csv(path_df, index = False)
    ####################################################################################################
    #now identify seq with the same cell_id (single cell data)... quite rare
    for _, df_cdr3_tmp in df_cdr3.iterrows():
        try:
            cdr3_TRA = df_cdr3_tmp['cdr3_TRA']
            cdr3_TRB = df_cdr3_tmp['cdr3_TRB']
            df_TRA = d_chain_seq[('TRA', cdr3_TRA)]
            df_TRB = d_chain_seq[('TRB', cdr3_TRB)]
            #keep only single cell data (with cell_id)
            df_TRA = df_TRA.dropna(subset = ['cell_id'])
            df_TRB = df_TRB.dropna(subset = ['cell_id'])
            df_paired = pd.merge(df_TRA, df_TRB, on = 'cell_id')
            if len(df_paired) > 0:
                key = [('TRA', cdr3_TRA, 'TRB', cdr3_TRB)]
                d_chain_seq[key] = df_paired
        except:
            continue

    #make directory and save results


print("\n Finished at --> {0} \n ".format(datetime.datetime.now()))


