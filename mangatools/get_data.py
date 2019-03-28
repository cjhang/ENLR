#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import sys
import datetime
import argparse
import requests
from bs4 import BeautifulSoup as bs
import re
from multiprocessing import Process

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target', help="target different files, available target: all, maps, datacube, qa, ref", type=str, default='all')
# parser.add_argument('-v', '--verbose', help="show the detail downloading imformation")
parser.add_argument('-p', '--plate', help="specify the plates to download, useage: -p '7443, 7597'")
args = parser.parse_args()

data_dir = 'MPL-6/'

if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

def download_file(url, target_path, size=None, auth=None, filename=None):
    # size is code in string
    s_local = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=10)
    s_local.mount(url, adapter)
    s_local.auth = auth
    from tqdm import tqdm
    is_download = False
    if not filename:
        filename = url.split('/')[-1]
    if filename in os.listdir(target_path):
        if size:
            # print('online sie is {}, local size is {}'.format(size, os.path.getsize(target_path+filename)))
            if size == str(os.path.getsize(target_path+filename)):
                print('Already downloaded {}, skipped...'.format(filename))
                return 0
    max_try = 5
    times_try = 0
    while True:
        if times_try > max_try:
            with open('get_data.log', 'a') as f:
                f.write("{}: Error enconter when download {}\n".format(datetime.datetime.now().isoformat(), filename))
            print("Failed! Cannot download {}, write to log file.".format(filename))
            return 0
        try:
            with s_local.get(url, stream=True, auth=auth, timeout=60) as r:
                handle = open(target_path+filename, "wb")
                tqdm.write("Downloading %s ..." % filename)
                for chunk in tqdm(r.iter_content(chunk_size=1024)):
                    if chunk:  # filter out keep-alive new chunks
                        handle.write(chunk)
            return
        except KeyboardInterrupt:
            sys.exit(1)
        except:
            print("Error encounter when download {}, re-downloading...".format(filename))
        times_try = times_try + 1

def download_ifugdsgn(data_dir, url_plateifu, session=None, target="all"):
    r_ifudsgn = session.get(url_plateifu)
    # print(r2.status_code)
    soup_ifudsgn = bs(r_ifudsgn.text, 'lxml').body.table.tbody
    links = [link.text for link in soup_ifudsgn.find_all('a')]
    sizes = [item.select('td')[1].text.strip() for item in soup_ifudsgn.select('tr')]
    # list contains all the job for multiprocessing
    # download_jobs = []

    if target == "datacube" or target == "all":
        # download logcube
        try:
            tmp = links[3]
        except KeyboardInterrupt:
            sys.exit(1)
        except:
            # split the url to get the plate and ifudsgn name
            datacube_name = data_dir.split('/')[-3] +'-'+ data_dir.split('/')[-2]
            with open('get_data.log', 'a') as f:
                f.write("{}: Error enconter when download manga-{}-LOGCUBE-SPX-GAU-MILESHC.fits.gz\n".format(datetime.datetime.now().isoformat(), datacube_name))
            print("Error, file not found! Cannot download datacube of {}, write to log file.".format(datacube_name))
            return
        download_file(url_plateifu+links[3], data_dir, size=sizes[3], auth=session.auth)
        # download_jobs.append((url_plateifu+links[3], data_dir, sizes[3], session.auth, None))
    if target == "maps" or target == "all":
        # download maps    
        # download_file(url_plateifu+links[4], data_dir, size=sizes[4], auth=session.auth)
        try:
            tmp = links[4]
        except KeyboardInterrupt:
            sys.exit(1)
        except:
            # split the url to get the plate and ifudsgn name
            map_name = data_dir.split('/')[-3] +'-'+ data_dir.split('/')[-2]
            with open('get_data.log', 'a') as f:
                f.write("{}: Error enconter when download manga-{}-MAPS-SPX-GAU-MILESHC.fits.gz\n".format(datetime.datetime.now().isoformat(), map_name))
            print("Error, file not found! Cannot download datacube of {}, write to log file.".format(map_name))
            return 0
        download_file(url_plateifu+links[4], data_dir, size=sizes[4], auth=session.auth)
        # download_jobs.append((url_plateifu+links[4], data_dir, sizes[4], session.auth, None))
    if target == "qa" or target == "all":
        # download all the files in qa
        if not os.path.isdir(data_dir+'qa/'):
            os.mkdir(data_dir+'qa/')
        qa_r = requests.get(url_plateifu + links[1], auth=session.auth)
        qa_soup = bs(qa_r.text, 'lxml').body.table.tbody
        qa_links = [link.text for link in qa_soup.find_all('a')]
        qa_sizes = [item.select('td')[1].text.strip() for item in qa_soup.select('tr')]
        for i in range(1,len(qa_links)):
            download_file(url_plateifu+links[1]+qa_links[i], data_dir+'qa/', qa_sizes[i], auth=session.auth)
            # download_jobs.append((url_plateifu+links[1]+qa_links[i], data_dir+'qa/', qa_sizes[i], session.auth, None))
    if target == 'ref' or target== "all":
        # download all the files in ref
        if not os.path.isdir(data_dir+'ref/'):
            os.mkdir(data_dir+'ref/')
        ref_r = requests.get(url_plateifu + links[2], auth=session.auth)
        ref_soup = bs(ref_r.text, 'lxml').body.table.tbody
        ref_links = [link.text for link in ref_soup.find_all('a')]
        ref_sizes = [item.select('td')[1].text.strip() for item in ref_soup.select('tr')]
        for i in range(1,len(ref_links)):
            download_file(url_plateifu+links[2]+ref_links[i], data_dir+'ref/', ref_sizes[i], auth=session.auth)
            # download_jobs.append((url_plateifu+links[2]+ref_links[i], data_dir+'ref/', ref_sizes[i], session.auth, None))
    # for job in download_jobs:
        # p = Process(target=download_file, args=job)
        # p.start()
        # p.join()
        # print('start {}'.format(job))

if __name__ == '__main__':
    
    url = "https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-6/SPX-GAU-MILESHC/"
    s = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=10)
    s.mount(url, adapter)
    s.auth = ('sdss', '2.5-meters')
    r = s.get(url)
    soup = bs(r.text, 'lxml')
    plates_online = [plate.get_text() for plate in soup.body.table.tbody.find_all('a')][1:]
    local_dir = os.listdir(data_dir)
    if args.plate:
        arg_plates = [plate.strip() for plate in args.plate.split(',')]
        plates = []
        for item in arg_plates:
            plate = item+'/' # plates_online were format like '7443/'
            if plate not in plates_online:
                print('Warning: plate-{} not fond in server, skipped'.format(item))
            else:
                plates.append(plate)
    else:
        plates = plates_online
    for plate in plates:
        if re.match(r'^\d+', plate).group() not in local_dir:
            os.mkdir(data_dir + plate)
        # download_plate(data_dir+plate, url+plate, session=s)
        plate_dir = data_dir+plate
        plate_url = url+plate
        r_plate = s.get(plate_url)
        soup_plate = bs(r_plate.text, 'lxml')
        ifudsgns = [ifudsign.get_text() for ifudsign in soup_plate.body.table.tbody.find_all('a')][1:]
        # print(ifudsigns)
        local_ifu = os.listdir(plate_dir)
        for ifudsgn in ifudsgns:
            if re.match(r'^\d+', ifudsgn).group() not in local_ifu:
                os.mkdir(plate_dir + ifudsgn)
            if args.target:
                download_ifugdsgn(plate_dir+ifudsgn, plate_url+ifudsgn, session=s, target=args.target)
            else:
                raise ValueError


