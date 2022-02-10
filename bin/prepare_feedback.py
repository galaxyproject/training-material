#!/usr/bin/env python

import datetime
import pandas as pd
from pathlib import Path

new_tutorial_name = {
    '16S Microbial Analysis with Mothur': '16S Microbial Analysis with mothur (extended)',
    '16S Microbial Analysis with mothur': '16S Microbial Analysis with mothur (extended)',
    'Quality Control ': 'Quality Control',
    'RNA-Seq reads to counts': '1: RNA-Seq reads to counts',
    'RNA-seq counts to genes': '2: RNA-seq counts to genes',
    'RNA-seq genes to pathways': '3: RNA-seq genes to pathways',
    'Introduction to Genome Assembly': 'An Introduction to Genome Assembly',
    'EWAS data analysis of 450k data': 'Infinium Human Methylation BeadChip',
    'Creating a new tutorial - Defining the technical infrastructure': 'Tools, Data, and Workflows for tutorials',
    'Creating a new tutorial - Writing content in Markdown': 'Creating content in Markdown',
    'Running the Galaxy Training material website locally': 'Running the GTN website locally',
    'Visualizations: JavaScript plugins': 'JavaScript plugins',
    'Compute and analyze Essential Biodiversity Variables with PAMPA toolsuite': 'Compute and analyze biodiversity metrics with PAMPA toolsuite',
    'Ephemeris for Galaxy Tool Management': 'Galaxy Tool Management with Ephemeris',
    'Collections: Rule Based Uploader': 'Rule Based Uploader',
    'Collections: Using dataset collection': 'Using dataset collections',
    'Data: Downloading and Deleting Data in Galaxy': 'Downloading and Deleting Data in Galaxy',
    'Histories: Understanding Galaxy history system': 'Understanding Galaxy history system',
    'Jupyter: Use Jupyter notebooks in Galaxy': 'Use Jupyter notebooks in Galaxy',
    'Using dataset collection': 'Using dataset collections',
    'Workflows: Extracting Workflows from Histories': 'Extracting Workflows from Histories',
    'Workflows: Using Workflow Parameters': 'Using Workflow Parameters',
    'Exome sequencing data analysis': 'Exome sequencing data analysis for diagnosing a genetic disease',
    'Galaxy Tool Management': 'Galaxy Tool Management with Ephemeris',
    'Virtual screening of the SARS-CoV-2 main protease with rDock and pose scoring': 'Virtual screening of the SARS-CoV-2 main protease with rxDock and pose scoring'
}
new_topic_for_tuto = {
    'Formation of the Super-Structures on the Inactive X': 'Epigenetics',
    'Identification of the binding sites of the Estrogen receptor': 'Epigenetics',
    'Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1)': 'Epigenetics',
    'RAD-Seq Reference-based data analysis': 'Ecology',
    'RAD-Seq de-novo data analysis': 'Ecology',
    'RAD-Seq to construct genetic maps': 'Ecology',
    'Advanced R in Galaxy': 'Foundations of Data Science',
    'R basics in Galaxy': 'Foundations of Data Science'
}
new_topics = {
    'User Interface and Features': 'Using Galaxy and Managing your Data',
    'Data Manipulation': 'Using Galaxy and Managing your Data',
    'User Interface and Data Manipulation': 'Using Galaxy and Managing your Data',
    'Assembly) is not working I can do up to multiQC and after unicycler not working': 'Assembly'
}

acceptable_topics = [
    "Assembly",
    "Climate",
    "Computational chemistry",
    "Contributing to the Galaxy Training Material",
    "Development in Galaxy",
    "Ecology",
    "Epigenetics",
    "Foundations of Data Science",
    "Galaxy Server administration",
    "Genome Annotation",
    "Imaging",
    "Introduction to Galaxy Analyses",
    "Metabolomics",
    "Metagenomics",
    "Proteomics",
    "Sequence analysis",
    "Statistics and machine learning",
    "Teaching and Hosting Galaxy training",
    "Transcriptomics",
    "Using Galaxy and Managing your Data",
    "Variant Analysis",
    "Visualisation"
]

def extract_tutorial_feedbacks(topic_df, topic_name):
    '''Extract pro/con per tutorial for a topic and
    write them in a file

    :topic_df: dataframe object for the topic
    :topic_name: name for the topic, name for the file
    '''
    grouped_by_tuto = topic_df.groupby(by="tutorial")
    with open('../results/%s.md' % topic_name, 'w') as f:
        for tuto, group in grouped_by_tuto:
            # get groups
            tuto_df = grouped_by_tuto.get_group(tuto)
            pros = []
            cons = []
            # get pros/cons
            for index, row in tuto_df.iterrows():
                if row['pro'] != 'nan':
                    pros.append("%s (*%s*)" % (row['pro'], row['timestamp']))
                if row['con'] != 'nan':
                    cons.append("%s (*%s*)" % (row['con'], row['timestamp']))
            # write in report file
            f.write("- **%s**\n" % tuto)
            if len(pros) > 0:
                f.write("  - Pro:\n    - ")
                f.write("\n    - ".join(pros))
            if len(cons) > 0:
                f.write("\n  - Con:\n    - ")
                f.write("\n    - ".join(cons))
            f.write("\n\n")


def fix_tutorial_info(df):
    '''Change tutorial topic or title

    :param df: dataframe
    '''
    df = df.copy()
    # change topic for some tutorials
    for tuto in new_topic_for_tuto:
        df.loc[df.tutorial == tuto, 'topic'] = new_topic_for_tuto[tuto]
    # rename topic for all tutorials in a topic
    for topic in new_topics:
        df.topic = (df
            .topic
            .replace(to_replace=topic, value=new_topics[topic]))
    # rename some tutorials
    for tuto in new_tutorial_name:
        df.loc[df.tutorial == tuto, 'tutorial'] = new_tutorial_name[tuto]
    return df


def extract_topic_tutorial_name(df):
    '''Extract topic from tutorial name

    :param df: dataframe
    '''
    df = df.copy()
    new = df['tutorial_topic'].str[::-1].str.split('(', n = 1, expand = True)
    df["tutorial"]= new[1].str[::-1].str[:-1]
    df["topic"]= new[0].str[::-1].str[:-1]
    return df


def prepare_feedback(url, out_file):
    '''Get and prepare feedback CSV file

    :param url: URL to Google sheet with feedback answers
    :param out_file: Path to output file
    '''
    df = (pd.read_csv(url, sep='\t')
        # rename column
        .rename(columns = {'Timestamp': 'timestamp',
                            'How much did you like this tutorial?': 'note',
                            'What did you like?': 'pro',
                            'What could be improved?': 'con',
                            'Tutorial': 'tutorial_topic',
                            'Your feedback is always anonymous. Also make it confidential (only visible to admins)?': 'anonymous'})
        # extract topic from tutorial name
        .pipe(extract_topic_tutorial_name)
        # remove rows with NaN on note, pro and con
        .dropna(subset=['note', 'pro', 'con'], how='all')
        # replace NaN in note by 0
        .fillna(value={'note': 0})
        # fill other NaN by empty string
        .fillna('')
        # format
        .assign(
            #note to integer
            note=lambda x: x['note'].astype(int),
            # format pro and con to string
            pro=lambda x: x['pro'].astype(str),
            con=lambda x: x['con'].astype(str),
            # extract month and date
            timestamp=lambda x: pd.to_datetime(x['timestamp'], dayfirst=True))
        .assign(
            month=lambda x: x['timestamp'].dt.strftime("%Y-%m"),
            date=lambda x: x['timestamp'].dt.strftime("%Y-%m-%d"))
        # change topic for some tutorials
        .pipe(fix_tutorial_info)
        # remove extra columns
        .drop(columns=["tutorial_topic", "timestamp"])
        # remove some rows
        .query('topic != "TESTIN"')
        .query('topic != "tes"')
        .query('topic != ""')
        # remove tutorials
        .query('tutorial != "RNA-seq counts to genes and pathways"')
        .query('tutorial != "Recording Job Metrics"')
        # Remove non-existent topics
        .query("topic in @acceptable_topics")
        #
        .query('pro != "TESTING"')
        .to_csv(out_file))


if __name__ == '__main__':
    url = 'https://docs.google.com/spreadsheets/d/1NfZhi5Jav7kl9zFCkeb7rIC2F8xW1isruv1TeO4WpNI/export?format=tsv'
    out_file = Path('metadata') / Path('feedback.csv')
    prepare_feedback(url, out_file)
