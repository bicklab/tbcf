#!/usr/bin/env python
# coding: utf-8

# In[9]:


import pandas as pd
#import polars as pl
import os
import glob
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess


# In[10]:


def download_files(file_urls, folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    os.chdir(folder_name)
    for url in file_urls:
        filename = url.split('/')[-1].strip()
        subprocess.run(['wget', '-O', filename, url], check=True)
    os.chdir('..')
    return os.path.abspath(folder_name)


# In[11]:


tcf_urls = ["https://dl.ew2.dnanex.us/F/D/YbJq1f2Yx2P6VZqQpJ23G3vXXx8pZ9ZJf57kqxp3/TCF_ukb23159_c2_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/BBb5q4G6YY284Y5yb5bZJQJZkgPy2Y0x1XBfq9P2/TCF_ukb23159_c1_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/y66Z5z570Jp8FPXkk7GXqFYX4bxf7pgqJY8q49y6/TCF_ukb23159_c6_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/kQ8V4XK26v7zbYBb42FbxK2Y5f0J8656jppxGGfk/TCF_ukb23159_c7_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/b8bPxV90VkF0B9P5yvfZjpKjVZV1XP7F2fxjVXjb/TCF_ukb23159_c8_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/xq397p47XGJvBVKX4z8xQ39FbZ5V9vFk8P7QxyP4/TCF_ukb23159_c18_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/g3zGpYX7BQvpFP0qvzJXzv6qyQBBVZGb3q2bgpx0/TCF_ukb23159_c3_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/ZYG1Q3fvfG0VGP78vVx6fz6kyKFBb355JjFXVZXg/TCF_ukb23159_c11_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/8zpv891KpVZK75JY7KY5KP45qF41vzK53F9xfQ7b/TCF_ukb23159_c19_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/JQ931X83pz5p4Gz19Bq6g1Vkk5b0y876G8jb468p/TCF_ukb23159_c17_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/jF4Bb2KfKbqKvkfG924kz5p2F1vyzF13pJjVvgYF/TCF_ukb23159_c10_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/f5j92y1Vz37ZJGp418yGb6Pz5bv1YfQ9f3X34yYP/TCF_ukb23159_c9_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/P2Xg3qzJ6xb4pBB5Z71q0YxXPyQGp8Q7Qxz3g3B4/TCF_ukb23159_c20_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/Vb6yK6QQ419jQ3kPyKGfBvVyqKfvGJGXfYPY8gVy/TCF_ukb23159_c5_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/64YGFQ5Jqfv6Z79VpKQv6p63ZKqP4kjz97178F1x/TCF_ukb23159_c15_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/Yp6fzyfFJVKF1BfF3pYGbQ1bv8ByqFx4fbfYg5Z8/TCF_ukb23159_c14_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/Yxx5BK73gQ4376qgGZ5zBPKVzkFy3x3pQPZYG744/TCF_ukb23159_c16_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/5z9Yy11Kq8Yy1yZ2z4680zK9KKx5QKJqVYgk88x7/TCF_ukb23159_c4_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/2xj32pqFykK2FV830JJq36z8kx603bGgxPKJF7FB/TCF_ukb23159_c12_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/qq3ZqgqPVx1ByP4QGvYbZj3ZX1fq0G699GFyPZgq/TCF_ukb23159_c13_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/qp591G0jz58P81kKxQBvQ9VQQQb6b0xKYk95gQY4/TCF_ukb23159_c22_b0_v1.bgen_TCF.regenie",
"https://dl.ew2.dnanex.us/F/D/Q5Yv4YBYg5Fbzg9K5FGB4BFjBGb5qvF468g573bQ/TCF_ukb23159_c21_b0_v1.bgen_TCF.regenie"]


# In[12]:


bcf_urls = ["https://dl.ew2.dnanex.us/F/D/x7Q391k6B8gX3KzZf8KfzjxpYY7P5F8pjjBFpJk7/BCF_ukb23159_c19_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/GFfVg3f66fqK8q5GVv4B1qV0vKvfp81y92k30V9x/BCF_ukb23159_c2_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/yFbf4F4Xq7KqZP36KbG4pFyPpzjqp3Vf3kPY223g/BCF_ukb23159_c11_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/xxJQk8pVKQ3GGp7y5J00YfQ0fqqz74ZXg509pB2f/BCF_ukb23159_c3_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/FpB6jf9Q5x34096zV0Zk1qjqPzfqzFxjJJ9K1VJ5/BCF_ukb23159_c6_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/598YQ158YX2jfqVy5Jbk5V4Zx14qz7xBQvK1BYFv/BCF_ukb23159_c10_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/bjP40zqkKpJF2vzzVBG2ByJV7gxFy95ZqY9PFK7z/BCF_ukb23159_c4_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/3v8P5pKVJV09G868KB0KjBvZP7722XF2yzK6K9qV/BCF_ukb23159_c17_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/f3gp6xyv1G62BVkpxvP05z86gv4pP0bq7k9bgZ44/BCF_ukb23159_c1_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/ZGGv4qZKgv8gKjx5KKZQx1B93z83vZFXv3x40FpK/BCF_ukb23159_c7_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/F04PzQGF7Fq0JVBggxFvq7j38K263xZzyXG50F1Z/BCF_ukb23159_c12_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/xpxyF2yp2y8K2jYjXzKfXYJbXb2ZqP33bPX1Bzyk/BCF_ukb23159_c16_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/yGBjzxG3Fxj4jPJvG14G8B2fGPp5V23KV3yyFqpJ/BCF_ukb23159_c15_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/qyp4FJvYYyXf6jqGV4GxPB97BGj8k420f0Pxg9x9/BCF_ukb23159_c9_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/yXGK1QYjGF88jG9Z71Kv11pfqvjVkZ8kjQvVvf6j/BCF_ukb23159_c14_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/3bpJP9kXQp2gy0J7jB91VJ2fq874ZVYXPy0GJJk2/BCF_ukb23159_c8_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/7ggZf2y5B8qkB1qKP3XVKFY0jfzKq48yy4Yk8GV0/BCF_ukb23159_c5_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/K1B7bJ440X3BJy2xgy88Pb2Zj9PYbB616Kb2X8yy/BCF_ukb23159_c20_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/1Gq1V1XYbfzz1jbFGgqXb27VF94x992X204k6b8P/BCF_ukb23159_c22_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/BGVpGKvg2qZV71YQ1VYxJ0BJbyZ6ZZfkkZ4093GQ/BCF_ukb23159_c18_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/FK5kg72ZBkkpzF5qfj34Kqpk5zJvvyv3KkGpfyyV/BCF_ukb23159_c13_b0_v1.bgen_BCF.regenie",
"https://dl.ew2.dnanex.us/F/D/KbZzqbB3J0J9JjPB95jKgQKYfF7JQg6x6GBZ2YVy/BCF_ukb23159_c21_b0_v1.bgen_BCF.regenie"]


# In[13]:


download_files(tcf_urls, "ukb_rare_tcf")


# In[14]:


download_files(bcf_urls, "ukb_rare_bcf")


# In[21]:


bcf_files =[]
for i in range(1,23):
    path = "ukb_rare_bcf/BCF_ukb23159_c" + str(i) + "_b0_v1.bgen_BCF.regenie"
    bcf = pd.read_table(path, skiprows=1, sep='\s+')
    bcf_fil = bcf[bcf['TEST']=='ADD-SKATO']
    bcf_files.append(bcf_fil)


# In[22]:


bcf_df = pd.concat(bcf_files, ignore_index=True)
bcf_df


# In[26]:


bcf_df.to_csv("gs://bicklab-main-storage/Users/Hannah_Poisner/ukb_bcf_rare_clean.regenie",index=False,sep='\t')


# In[27]:


tcf_files =[]
for i in range(1,23):
    path = "ukb_rare_tcf/TCF_ukb23159_c" + str(i) + "_b0_v1.bgen_TCF.regenie"
    tcf = pd.read_table(path, skiprows=1, sep='\s+')
    tcf_fil = tcf[tcf['TEST']=='ADD-SKATO']
    tcf_files.append(tcf_fil)


# In[28]:


tcf_df = pd.concat(tcf_files, ignore_index=True)
tcf_df


# In[29]:


tcf_df.to_csv("gs://bicklab-main-storage/Users/Hannah_Poisner/ukb_tcf_rare_clean.regenie",index=False,sep='\t')


# In[ ]:




