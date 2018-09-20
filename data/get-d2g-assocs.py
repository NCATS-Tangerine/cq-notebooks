import requests


# Globals and Constants
SOLR_URL = 'https://solr-dev.monarchinitiative.org/solr/golr/select'


def main():

    result_docs = get_causal_disease_gene_assocs()
    get_disease_phenotype_list(result_docs)

    for key, value in result_docs.items():
        print("{}\t{}\t{}\t{}\t{}".format(
            key,
            value['label'],
            len(value['gene_ids']),
            value['phenotype_count'],
            "|".join(value['gene_labels'])
        ))


def get_disease_phenotype_list(result_docs):

    for key, val in result_docs.items():

        filters = ['subject_closure:"{0}"'.format(key),
                   'object_category:"phenotype"'
                   ]
        params = {
            'wt': 'json',
            'rows': 0,
            'start': 0,
            'q': '*:*',
            'fq': filters,
            'facet': 'true',
            'facet.mincount': 1,
            'facet.sort': 'count',
            'json.nl': 'arrarr',
            'facet.limit': -1,
            'facet.field': 'object'
        }
        solr_request = requests.get(SOLR_URL, params=params)
        response = solr_request.json()
        result_docs[key]['phenotype_count'] = len({val[0] for val in response['facet_counts']['facet_fields']['object']})

    return result_docs


def get_causal_disease_gene_assocs():
    print("Fetching causal human gene phenotype and disease associations")
    result_docs = {}
    filters = ['object_category:"{0}"'.format("disease"),
               'subject_category:"gene"',
               'subject_taxon: "{0}"'.format('NCBITaxon:9606')]
    params = {
        'wt': 'json',
        'rows': 1000,
        'start': 0,
        'q': '*:*',
        'fq': filters,
        'fl': 'subject, object, object_label, subject_label, relation, is_defined_by'
    }
    
    causal_source = ["https://data.monarchinitiative.org/ttl/clinvar.ttl",
                     "https://data.monarchinitiative.org/ttl/omim.ttl",
                     "https://data.monarchinitiative.org/ttl/orphanet.ttl"]
    resultCount = params['rows']
    while params['start'] < resultCount:
        solr_request = requests.get(SOLR_URL, params=params)
        response = solr_request.json()
        resultCount = response['response']['numFound']

        for doc in response['response']['docs']:
            if 'relation' in doc:
                # Filter out likely pathogenic
                if doc['relation'] == 'GENO:0000841':
                    continue

            if 'is_defined_by' in doc\
                    and len([source for source in doc['is_defined_by'] if source in causal_source]) == 0\
                    and doc['is_defined_by'] != ['https://data.monarchinitiative.org/ttl/hpoa.ttl']:
                    continue

            if doc['object'] in result_docs:
                result_docs[doc['object']]['gene_ids'].append(doc['subject'])
                result_docs[doc['object']]['gene_labels'].append(doc['subject_label'])
            else:
                result_docs[doc['object']] = {
                    'label': doc['object_label'],
                    'gene_ids': [doc['subject']],
                    'gene_labels': [doc['subject_label']],
                    'phenotype_count': 0
                }

        params['start'] += params['rows']

    return result_docs


if __name__ == "__main__":
    main()
