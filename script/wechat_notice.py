#!/home/donghj/snakemake/bin/python
#wechat_notice.py
#version: 0.2
#pipeline: none
#author: dhj
#update: 2025/07/31

import argparse
import requests
import json
import os

#text color : comment:grey info:green warning:orange
upload_url= r'https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=d2b0093e-2be4-4fe2-a379-543a3116d37c'

def message(text,text_color):
    tagged_text="\n".join(["<font color=\""+text_color+"\">"+txt+"</font>" for txt in text.split("\n")])
    return(tagged_text)
def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--status',required=True,type=str, help='message status')
    parser.add_argument('--AT',nargs='+',type=str, default=['所有人'], help='AT people on wechat')
    parser.add_argument('--task',default='no_id',type=str, help='task number of this message')
    parser.add_argument('--log',default='',type=str, help='log file of message')
    args = parser.parse_args()
    
    location=os.getcwd()
    if args.status == "success":
        text=message("Congratulations! Task ",'info')+args.task+message(" successfully finished.\n",'info')
    elif args.status == "error":
        text=message("Oops! There're some trouble with task ",'warning')+args.task+message(", please check it.\n",'warning')
    else:
        text=message("Just a normal notification about task ",'comment')+args.task+message(", consider using \"success\" or \"error\" for more explicit reminder.\n",'comment')
    upload_text={
        "msgtype": "markdown",
        "markdown": {
            "content": " ".join(["<@"+name+">" for name in args.AT])+"\n"+text+message(args.log,'comment')+"\nProject Path: "+location
        }
    }

    response = requests.post(upload_url, headers={'Content-Type': 'application/json'}, data= json.dumps(upload_text))

if __name__ == '__main__':
    main() 





