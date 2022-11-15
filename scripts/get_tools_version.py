#!/usr/bin/env python

def get_dependencies(file) :
    """Return a list of the tools and versions from an env.yaml file."""
    tools = []
    dep = False
    with open(file, 'r') as f :
        for line in f :
            if dep :
                tools.append(line.replace(" ", "")[1:-1])
            if line == "dependencies:\n" :
                dep = True
    return tools


def get_CLR_tools(env_path) :
    """Return list of tools and version used for CLR data."""
    tools = []
    pbsv = env_path + "/pbsv_env.yaml"
    tools += get_dependencies(pbsv)
    longshot = env_path + "/longshot_env.yaml"
    tools += get_dependencies(longshot)
    tabix = env_path + "/tabix_env.yaml"
    tools += get_dependencies(tabix)
    return tools


def get_hifi_tools(env_path) :
    """Return list of tools and version used for hifi data."""
    tools = []
    pbsv = env_path + "/pbsv_env.yaml"
    tools += get_dependencies(pbsv)
    tools += ["deepvariant:singlularity-3.5.3"]
    tabix = env_path + "/tabix_env.yaml"
    tools += get_dependencies(tabix)
    return tools


def get_ont_tools(env_path) :
    """Return list of tools and version used for ont data."""
    tools = []
    svim = env_path + "/svim_env.yaml"
    tools += get_dependencies(svim)
    tools += ["pepper:singlularity-3.5.3"]
    tabix = env_path + "/tabix_env.yaml"
    tools += get_dependencies(tabix)
    return tools


def get_tools(env_path, data_type) :
    """Return list of tools and version used depending on the data."""
    if data_type == 'CLR' :
        return get_CLR_tools(env_path)
    if data_type == 'CCS' or data_type == 'hifi' :
        return get_hifi_tools(env_path)
    if data_type == 'ONT' :
        return get_ont_tools(env_path)


# print(get_tools("/home/quentinboone/Bureau/v16/envs", 'CLR'))
