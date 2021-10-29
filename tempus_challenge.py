#! /usr/local/bin/python3

### load libraries
import subprocess
import os
import requests
import pandas as pd
import argparse


### Part 1: Exploratory analysis
### Understand the current contents on the input VCF
def part1_exploratory_analysis(vcf,prefix,verbose):
    if verbose == True:
        print('Performing exploratory analysis.\n')
    
    
    # What information is in the INFO column?
    if verbose == True:
        print('What information is in the INFO column?\n')
        
    command = " grep '##INFO' {0} | cut -f3,6 -d'=' | sed 's/,Number=/\t/g' | sed 's/>//g'".format(vcf)
    output = subprocess.check_output(command,shell=True).decode('UTF-8')
    if verbose == True:
        print(output + '\n')
    
    
    # What information is in the FORMAT column?
    if verbose == True:
        print('What information is in the FORMAT column?\n')
    command = "grep '##FORMAT' {0} | cut -f3,6 -d'=' | sed 's/,Number=/\t/g' | sed 's/>//g'".format(vcf)
    output = subprocess.check_output(command,shell=True).decode('UTF-8')
    if verbose == True:
        print(output + '\n')
        
        
    # How many variants are in the file?
    if verbose == True:
        print('How many variants are in the file?\n')
    variant_count = "grep -v '#' {0} | wc -l".format(vcf)
    variant_count = int(subprocess.check_output(variant_count,shell=True).decode('UTF-8'))
    if verbose == True:
        print('{0} variants are in the input VCF file.\n'.format(variant_count))
    
    # How many samples are in the file?
    if verbose == True:
        print('How many samples are in the file?\n')
    # command which uses bcftools to get sample IDs from the file
    command = "bcftools query -l {0}".format(vcf)
    # list of samples in vcf
    vcf_samples = subprocess.check_output(command,shell=True).decode('UTF-8').strip().split('\n')
    if verbose == True:
        print('These are the samples in the input VCF file: '+', '.join(vcf_samples) + '\n')
    
    return





### Read in the VCF
# Read the vcf into a pandas dataframe. Keep the comments saved into a separate list. 
# Additional comments generated will be added to the comment list, and in the end the
# comment list will be merged with the updated VCF body. 
def part2_read_VCF(vcf,prefix,verbose):
    
    if verbose == True:
        print('Reading the VCF\n')
    
    # command to get the number of lines which are comment only (comment lines begin with double hashes)
    # we need the number of comment lines in order to know how many lines to skip when reading in the VCF into a df
    comment_count = "grep '##' {0} | wc -l".format(vcf)
    comment_count = int(subprocess.check_output(comment_count,shell=True).decode('UTF-8'))
    if verbose == True:
        print('The number of comments is: {0}\n'.format(comment_count))

    # read the vcf into a pandas dataframe, skipping the comment lines
    vcf_df = pd.read_csv(vcf,sep='\t',skiprows=comment_count,header=0)
    
    
    # save the original comments to a list
    command = "grep '##' {0}".format(vcf)
    comments_list = subprocess.check_output(command,shell=True).decode('UTF-8').strip().split('\n')
    
    
    
    # Add variant ID to the VCF. 
    # This is not the dbSNP rsID, but for now it is a description of the locus. 
    # The format for the ID is CHR-POS-REF-ALT.
    vcf_df['ID'] = vcf_df['#CHROM'].astype(str) + '-' + vcf_df['POS'].astype(str) + '-' + vcf_df['REF'] + '-' + vcf_df['ALT']

    
    
    return vcf_df,comments_list





### Annotate VCF: Allele frequency, variant type, variant impact, most deleterious variant impact
# For multiple consequences, annotate with the most deleterious possibility.
# In the ExAC database, the most deleterious consequence comes first. 
def part3_annotate_VCF(vcf_df,comments_list,prefix,verbose):
    
    if verbose == True:
        print('Annotating variant\n')
        
    # initialize a list of INFO tags
    info_tags = []

    # iterate through the variant IDs
    for variant_ID in vcf_df['ID']:

        # API URL to get the variant json
        url = 'http://exac.hms.harvard.edu/rest/variant/variant/{0}'.format(variant_ID)

        # pull the information
        response = requests.get(url)

        # convert to json
        variant_info = response.json()

        # initialize the list of consequences, and major consequences, biotype
        consequences_set = set()
        major_consequences_set = []
        biotype_set = set()

        # complete information is not available for every variant (e.g. allele frequency), 
        # so extract information through try-except
        try: 
            # extract allele frequency, round it to 3 decimal spaces
            allele_freq = round(variant_info['allele_freq'],3)


            # extract vep annotations, which contain the variant consequence
            # vep_anno is a list of dictionaries, each dictionary specifies a different version of the ___
            vep_anno = variant_info['vep_annotations']

            # check to see if vep_anno is empty of annotations
            if vep_anno == []:

                # set variables to '.'
                consequences_set.add('.')
                major_consequences_set.append('.')
                biotype_set.add('.')

            else:

                # iterate through the versions in vep_anno
                # note: the most deleterious consequence will be in the first iteration in the list of vep_anno,
                # according to: https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
                for version in vep_anno:

                    # extract consequence(s)
                    consequences = version['Consequence']

                    # extract the major consequence 
                    ### THE first one in the iteration will be the major
                    major_consequence = version['major_consequence']

                    ### ADDITIONAL annotation ###
                    # extract the biotype
                    biotype = version['BIOTYPE']

                    # if any of the values are equal to the empty string, set them equal to '.' instead
                    if biotype == '':
                        biotype = '.'
                    if major_consequence == '':
                        major_consequence = '.'
                    if consequences == '':
                        consequences = '.'

                    # add consequence to list
                    consequences_set.add(consequences.replace('&',','))
                    # add major consequence to the set
                    major_consequences_set.append(major_consequence)
                    # add biotype to set
                    biotype_set.add(biotype)

        except:
            # set variables to '.'
            allele_freq = '.'
            consequences_set.add('.')
            major_consequences_set.append('.')
            biotype_set.add('.')


        # make a tag for frequence
        freq_tag = 'ExAC_AF={0}'.format(allele_freq)


        # convert all the consequences into a string delimited with a comma
        consequences_str = ','.join(set(','.join(consequences_set).split(',')))
        # make a tag for consequences
        consequence_tag = 'CSQ={0}'.format(consequences_str)


        # the major consequece is the first in the major consequence set
        major_consequence = major_consequences_set[0]
        # make a tag for major consequence
        major_consequence_tag = 'Major_CSQ={0}'.format(major_consequence)


        # convert set of biotypes to str delimited with a comma
        biotype_str = ','.join(biotype_set)
        # make a tag for biotype
        biotype_tag = 'BIOTYPE={0}'.format(biotype_str)


        # concatenate the tags
        tag = freq_tag + ';' + consequence_tag + ';' + major_consequence_tag + ';' + biotype_tag
        # add tag to list of tags
        info_tags.append(tag)
    
    
    # Add the additional INFO tags to the VCF data frame's INFO column
    tags = pd.DataFrame({'TAG':info_tags})
    vcf_df['INFO'] = vcf_df['INFO']+';'+ tags['TAG']
    
    
    # For the three tags, Exac_AF, CSQ, Major_CSQ, and BIOTYPE, add comment lines to the comment list. 
    # make a comment line to explain the tags
    exac_af_comment = '##INFO=<ID=ExAC_AF,Number=A,Type=Float,Description="The ExAC allele frequency of the variant.">'
    consequence_comment = '##INFO=<ID=CSQ,Number=A,Type=String,Description="The consequence(s) of the variant, such as missense variant, inframe deletion, etc.">'
    major_consequence_comment = '##INFO=<ID=Major_CSQ,Number=A,Type=String,Description="The most severe consequence, as defined by Ensembl (https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences).">'
    biotype_comment = '##INFO=<ID=BIOTYPE,Number=A,Type=String,Description="A gene or transcript classification, as defined by Ensembl (https://m.ensembl.org/info/genome/genebuild/biotypes.html).">'

    
    # add comments to list of comments
    comments_list.append(exac_af_comment)
    comments_list.append(consequence_comment)
    comments_list.append(major_consequence_comment)
    comments_list.append(biotype_comment)


    # Annotate VCF with percentage of reads supporting the variant
    # split the INFO line by tags (;), then separate the key and value (=), and extract the value [1]
    # total read depth
    dp = vcf_df['INFO'].str.split(';',expand=True)[7].str.split('=',expand=True)[1].astype(int)

    # total alt var read observation count
    # note: in some cases more than one number is provided, so that has to be parsed later
    ao = vcf_df['INFO'].str.split(';',expand=True)[5].str.split('=',expand=True)[1]

    # total ref var read observation count 
    ro = vcf_df['INFO'].str.split(';',expand=True)[28].str.split('=',expand=True)[1]


    # make a sum of the alt var observation counts when more than one number is provided
    # initialize the list of alt var read observation counts
    ao_summed_list = []

    # iterate through each row
    for ao_row in ao.str.split(','):
        # converts the value into a list of integers
        ao_row = [int(ao_value) for ao_value in ao_row]

        # sums the integers (important when more than one integer)
        ao_row = sum(ao_row)

        # adds to list
        ao_summed_list.append(ao_row)

    # make a dataframe of the alt var observation counts
    ao_df = pd.DataFrame({'ID':vcf_df['ID'],'AO':ao_summed_list})

    # calculate the percent reads which support the variant
    # AO / (AO + RO)
    percentage = round(ao_df['AO'] / (ao_df['AO'] + ro.astype(int)) * 100,2)

    # create tag for percent reads which support the variant
    tag = 'PSV=' + percentage.astype(str)

    # create comment for percent reads which support the variant
    comment_psv = '##INFO=<ID=PSV,Number=1,Type=Float,Description="Percentage of reads supporting the variant versus those supporting reference reads.">'
    # add comment to comment list
    comments_list.append(comment_psv)

    # add PSV tag to the VCF INFO column
    vcf_df['INFO'] = vcf_df['INFO'] + ';' + tag
    
    
    # Write the VCF dataframe to a text file, and the updated comments to another text file, and merge the text files.
    if verbose == True:
        print('Writing annotated VCF to file\n')
            
    # vcf body file name
    vcf_body = prefix + '_body.txt'
    # write to a file
    vcf_df.to_csv(vcf_body,index=False,sep='\t')

    # vcf comments file name
    vcf_comments = prefix + '_comments.txt'
    # convert comments list to a string
    comments_string = '\n'.join(comments_list) + '\n'
    # write to a file
    handle = open(vcf_comments,'w')
    handle.write(comments_string)
    handle.close()

    # merge the comments file and the VCF body together. 
    vcf_updated = prefix + '_annotated.vcf'
    command = "cat {0} {1} > {2} ".format(vcf_comments,vcf_body,vcf_updated)
    subprocess.call(command,shell=True)

    # remove the body and comments file
    command = "rm {0} {1}".format(vcf_comments, vcf_body)
    subprocess.call(command,shell=True)
    
    return 
    
    
    



### Main function
if __name__ == "__main__":

    
    ### set up command for taking in command line arguments
    parser = argparse.ArgumentParser()


    ### define what arguments to take
    ### VCF input and output prefix are required 
    # VCF file
    parser.add_argument('-i','--input', type=str,help="input VCF file",required=True)
    # output file prefix
    parser.add_argument('-p','--prefix', type=str,help="output VCF file prefix",required=True)
    # verbose mode
    parser.add_argument("-v", "--verbose", action="store_true",help="turn verbose mode off",default=True)
    
    
    # parses arguments
    args = parser.parse_args()
    
    
    # extract arguments
    vcf = args.input
    prefix = args.prefix
    verbose_mode = args.verbose
    
    # call functions
    # expoloratory analysis
    part1_exploratory_analysis(vcf,prefix,verbose_mode)
    # read VCF
    vcf_df,comments_list = part2_read_VCF(vcf,prefix,verbose_mode)
    # annotate VCF
    part3_annotate_VCF(vcf_df,comments_list,prefix,verbose_mode)
