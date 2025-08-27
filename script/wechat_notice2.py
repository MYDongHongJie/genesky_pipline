#!/home/genesky/software/python/3.9.4/bin/python3
#_*_coding:utf-8 _*_

import argparse
import datetime
import getpass
import glob
import json
import os
import re
import socket
import ssl
import sys
import urllib
from urllib import request

import simplejson

# def read_quality_check_file(*files):
#     """读入并合并所有的指控指标文件。每个文件的表头都要遵循 sample  指标1  指标2  ...  QC 这种格式
#     Returns:
#         _type_: 返回两个数据，1. 字典，记录每个样本的详细指标，以及最终汇总的QC. 2. 指标列表
#     """
#     sample_quality = {}
#     indicators = []  # 检测指标名称
#     for file in files:
#         with open(file, 'r') as fh:
#             titles = fh.readline().strip().split('\t')
#             indicators.extend(titles[1:-1])
#             for line in fh:
#                 values = line.strip().split('\t')
#                 if line.startswith('#') or len(values) != len(titles):
#                     continue
#                 sample, qc = values[0], values[-1]
#                 # 初始化
#                 if sample not in sample_quality:
#                     sample_quality[sample] = {'QC_merge': []}
#                 # 每一个质控文件的qc信息放在一起, PASS 先不放
#                 for qc_value in qc.split(','):
#                     if qc_value != 'PASS':
#                         sample_quality[sample]['QC_merge'].append(qc_value)
#                 # 数据记录
#                 for title, value in zip(titles, values):
#                     sample_quality[sample][title] = value
#     # 填充缺失值/QC汇总
#     total_qc = ''
#     for sample in sample_quality.keys():
#         for indicator in indicators:
#             if indicator not in sample_quality[sample]:
#                 sample_quality[sample][indicator] = '.'
#         # qc
#         is_error = len([qc for qc in sample_quality[sample]['QC_merge'] if 'ERROR' in qc]) > 0
#         if not is_error:
#             sample_quality[sample]['QC_merge'].insert(0, 'PASS')
#         sample_quality[sample]['QC'] = ",".join(sample_quality[sample]['QC_merge'])
#         total_qc += sample_quality[sample]['QC']

#     condition = '-1'  # -1 表示未知
#     if 'PASS' in total_qc:
#         condition = '0'
#     if 'WARNING' in total_qc:
#         condition = '2'
#     if 'ERROR' in total_qc:
#         condition = '1'
    
#     # 整理成文本
#     text_all = "sample\t" + "\t".join(indicators) + "\tQC\n"
#     text_error = "sample\t" + "\t".join(indicators) + "\tQC\n"
#     text_warning = ""
#     for sample in sample_quality.keys():
#         values = [sample] + [ sample_quality[sample][indicator] for indicator in indicators] + [sample_quality[sample]['QC']]
#         text_all += "\t".join(values) + "\n"
#         if 'ERROR' in sample_quality[sample]['QC']:
#             text_error += "\t".join(values) + "\n"
#         elif 'WARNING' in sample_quality[sample]['QC']:
#             text_warning += "\t".join(values) + "\n"

#     return sample_quality, indicators, text_all, text_error, text_warning, condition


# def check_is_new_form(files):
#     is_new_form = True
#     for file in files:
#         with open(file, 'r') as fh:
#             titles = fh.readline().strip().split('\t')
#             # 表头最后一列必须是OK字符
#             if titles[-1] != 'QC':
#                 is_new_form = False
#             # 内容带有condition字符，也是旧版本
#             for line in fh:
#                 if 'condition' in line.lower():
#                     is_new_form = False
#     return is_new_form

def read_notify_qc_file(files):
    content = ''
    for file in files:
        with open(file, 'r') as fh:
            content += fh.read()
    return content

def GetOptions(args = sys.argv[1:]):
    parser = argparse.ArgumentParser(description = "根据项目路径下的文件推算出项目,用于项目完成之后发送微信提醒")
    parser.add_argument("-i", "--qc_files", required=True, help = "qc_files")
    parser.add_argument("-c", "--project", help ="项目号")
    parser.add_argument("--platforms", help ="10X or huadaC4")
    parser.add_argument("--result_dir", help ="结果路径")
    options = parser.parse_args(args)
    return options

def gettoken(corpid,corpsecret):
    context = ssl._create_unverified_context()
    gettoken_url = 'https://qyapi.weixin.qq.com/cgi-bin/gettoken?corpid=' + corpid + '&corpsecret=' + corpsecret
    #  print  gettoken_url
    try:
        token_file = request.urlopen(gettoken_url, context=context)
    except urllib.error.HTTPError as e:
        print(e.code)
        print(e.read().decode("utf8"))
        sys.exit()
    token_data = token_file.read().decode('utf-8')
    token_json = json.loads(token_data)
    token_json.keys()
    token = token_json['access_token']
    return token

def senddata(access_token, send_username, message):
    context = ssl._create_unverified_context()
    send_url = 'https://qyapi.weixin.qq.com/cgi-bin/message/send?access_token=' + access_token
    send_values = {
        "touser":send_username,       #企业号中的用户帐号.
        # "toparty":"25",             #企业号中的部门id.
        "msgtype":"text",
        "agentid":"1000002",       #企业号中的应用id.
        "text":{
            "content":message
           },
        "safe":"0"
        }
    send_data = simplejson.dumps(send_values, ensure_ascii=False).encode('utf-8')
    # print(send_data)
    send_request = request.Request(send_url, send_data)
    response = json.loads(request.urlopen(send_request, context=context).read())
    return response

# 解析qc文件
def parse_qc_txt(txt):
    qc_condition = ''
    qc_detail = ''
    with open(txt, 'r') as fh:
        for line in fh:
            line = line.strip()
            if re.search(r'\w', line):
                if(re.search('condition', line, re.I)):
                    name, qc_condition = line.replace(' ', '').split(':', 1)
                else:
                    qc_detail += line + '\r\n'

    return qc_condition, qc_detail

def get_project_qc_info(qc_files):
    qc_info = {}
    qc_info['QC_DETAIL'] = ''
    qc_info['QC_CONDITION'] = 'QC UNKNOWN'
    # 解析qc.txt
    condition_list = []
    qc_condition, qc_detail = parse_qc_txt(qc_files)
    condition_list.append(qc_condition)
    qc_info['QC_DETAIL'] += qc_detail
    if '1' in condition_list:
        qc_info['QC_CONDITION'] = 'QC ERROR'
    elif '2' in condition_list:
        qc_info['QC_CONDITION'] = 'QC WARNING'
    elif '0' in condition_list:
        qc_info['QC_CONDITION'] = 'QC PASS'
    return qc_info

def split_message(message, limit=2048):
    '''一条消息的长度不能超过 limit 个字节'''
    messages = []
    tmp = ''
    length = 0
    for line in message.split('\n'):
        if (len(line) + length + 1) < limit:
            tmp += line + "\n"
        else:
            messages.append(tmp)
            tmp = line
        length = len(tmp)
    if tmp:
        messages.append(tmp)

    # 清理，防止空白
    messages_clean = []
    for m in messages:
        m = m.strip()
        if m:
            messages_clean.append(m)
    return messages_clean

if __name__ == '__main__':
    #企业号的标识ID.
    corpid     = 'wwb4584b8032301692'
    #应用程序的密钥
    corpsecret = 'zHkaAoD9aZG25Jz4RTLN5COcV7PXDvZdrAB0J2926bw'
    options = GetOptions()

    user_id = {'lch':'lch','ganb':'ganb','xuy':'xux', 'zhouyh':'zhouyh', 'wangr': 'wangr', 'guansb': 'guansb', 'xiewt': 'xiewt', 'hanxh': 'hanxh', 'lingp': 'lingp', 'duzh': 'duzhihao', 'hexy': 'hexinyu', 'yangd': 'yangdan', 'lijk': 'lijiake', 'panyy': 'panyingying', 'donghj': 'donghongjie'}

    linux_user = getpass.getuser()
    user_name = "|".join([user_id[linux_user], user_id['lch'], user_id['ganb'],  user_id['guansb'] ])
    server_ip = socket.gethostbyname(socket.gethostname())

    # 获取项目路径、项目名称
    qc_files = str(options.qc_files)
    project = str(options.project)
    big_direction =  str(options.platforms)
    result_dir = str(options.result_dir)
    # 消息接收的用户新增：微生物方向的发给杨丹，其他发给周元昊
    if big_direction == 'MGS':
        user_name += "|" + user_id['panyy']
    elif big_direction == 'WTS' or big_direction == '10X' or big_direction == 'huadaC4':
        user_name += "|" + user_id['zhouyh']
    else:
        user_name += "|" + user_id['yangd']

    
    # 获取项目QC情况
    qc_info = get_project_qc_info(qc_files)
   
    # 获取企业微信访问密钥
    accesstoken = gettoken(corpid,corpsecret)
    # 发送消息 
    message_raw = f"{project}    {server_ip}    {linux_user}    finished ({qc_info['QC_CONDITION']})\n {result_dir} \n\n{qc_info['QC_DETAIL']}"
    # 拆分消息，分批发送，单条长度不超过 2048 （企业微信限制）
    messages = split_message(message_raw, 2048)
    for message in messages:
        response    = senddata(accesstoken, user_name, message)

        # 如果发送失败，则用邮件发送
        if response['errcode'] != 0:
            email    = "%s@geneskies.com" % linux_user
            Subject  = "[ERROR] %s send message failed!!!" % sys.argv[0]
            os.system("echo \"%s\" | mail -s \"%s\" %s" % (message, Subject, email))

