import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import seaborn as sns
from statsmodels.stats import multitest
from scipy import stats
from intervaltree import IntervalTree, Interval
from khipu.epdsConstructor import epdsConstructor
from jms.io import read_table_to_peaks
from jms.dbStructures import ExperimentalEcpdDatabase, knownCompoundDatabase
from khipu.extended import *
from khipu.utils import *
import shutil
import statsmodels.api as sm
from statsmodels.formula.api import ols

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerBase
from pca import pca
import sys



# read input data
feature_table = pd.read_csv(sys.argv[1], index_col=0, sep="\t")
table_name = sys.argv[2]
metadata = pd.read_csv(sys.argv[3], sep=",")
empCpds = json.load(open(sys.argv[4]))
annot_source = json.load(open(sys.argv[5]))
interactive = False
save_figs = True

formula_name_map = {}
formula_category_map = {}
formula_class_map = {}
category_class_annots = False
for cpd in annot_source:
    if cpd['neutral_formula'] not in formula_name_map:
        formula_name_map[cpd['neutral_formula']] = set()
        formula_class_map[cpd['neutral_formula']] = set()
    if cpd['name']:
        formula_name_map[cpd['neutral_formula']].add(cpd['name'])
    if 'category' in cpd and cpd['category']:
        if cpd['neutral_formula'] not in formula_category_map:
            formula_category_map[cpd['neutral_formula']] = set()
        category_class_annots = True
        formula_category_map[cpd['neutral_formula']].add(cpd['category'])
    if 'class' in cpd and cpd['class']:
        if cpd['neutral_formula'] not in formula_class_map:
            formula_class_map[cpd['neutral_formula']] = set()
        formula_class_map[cpd['neutral_formula']].add(cpd['class'])


feature_id_to_emp_cpd_id_map = {}
feature_id_to_cpd_name = {}
feature_id_to_category = {}
feature_id_to_class = {}
feature_id_to_emp_id = {}
for emp_id, emp_cpd in empCpds.items():
    for peak in emp_cpd["MS1_pseudo_Spectra"]:
        formulas = []
        names = []
        categories = set()
        classes = set()
        if 'list_matches' in emp_cpd:
            for match in emp_cpd['list_matches']:
                formula = match[0].split("_")[0]
                names_for_formula = formula_name_map[formula]
                formulas.append(formula)
                names.append(','.join(names_for_formula))
                if category_class_annots:
                    categories = categories.union(formula_category_map[formula])
                    classes = classes.union(formula_class_map[formula])
            feature_id_to_emp_cpd_id_map[peak['id_number']] = ','.join(formulas)
            feature_id_to_cpd_name[peak['id_number']] = ','.join(names)
            feature_id_to_category[peak['id_number']] = categories
            feature_id_to_class[peak['id_number']] = classes
            feature_id_to_emp_id[peak['id_number']] = emp_id
        #feature_id_to_emp_cpd_id_map[peak['id_number']] = ','.join([x[0].split("_")[0] for x in emp_cpd["list_matches"]]) if "list_matches" in emp_cpd else emp_cpd["neutral_formula_mass"]
        #feature_id_to_cpd_name[peak['id_number']] = ','.join([','.join(formula_name_map[x[0].split("_")[0]]) for x in emp_cpd["list_matches"]]) if "list_matches" in emp_cpd else emp_cpd["neutral_formula_mass"]

feature_table['annot'] = feature_table['id_number'].map(feature_id_to_emp_cpd_id_map)
feature_table['annot_names'] = feature_table['id_number'].map(feature_id_to_cpd_name)
feature_table['emp_id'] = feature_table['id_number'].map(feature_id_to_emp_id)
print(feature_table['emp_id'].values)
if category_class_annots:
    feature_table['categories'] = feature_table['id_number'].map(feature_id_to_category)
    feature_table['classes'] = feature_table['id_number'].map(feature_id_to_class)


for x in sys.argv[6:]:
    feature_table.drop(columns=x, inplace=True)

# find names of samples
all_sample_names = metadata["Name"]
sample_names = list(set(metadata["Name"]).intersection(feature_table.columns))

# fill in missing values
def calc_interpolate_value(row, indices):
    values = [x for x in list(row[indices]) if x > 0]
    if values:
        return np.min(values) / 2
    else:
        return 0
feature_table["feature_interpolate_value"] = feature_table.apply(calc_interpolate_value, axis=1, args=(sample_names,))
for sample_name in sample_names:
    feature_table[sample_name] = feature_table[[sample_name, "feature_interpolate_value"]].max(axis=1)


# set up PCA config

media_fill = {
    'basal': None,
    'complex': 'full'
}

differentation_color = {
    'ipsc': 'k',
    'differentiated': 'r'
}

genotype_marker = {
    'cyp19a1 KO': '*',
    'wt': 'o',
    'prlr KO': '^',
    'ptn KO': 'v'
}

additives_str = {
    "prolactin": '+p',
    "dheas": '+d',
    "None": ''
}

# do PCA
if interactive or save_figs:
    scaler = StandardScaler()
    transformed = scaler.fit_transform(np.log2(feature_table[sample_names].T+1))
    model = PCA(n_components=2)
    pca_transformed = model.fit_transform(transformed)
    for x,y,name in zip([x[0] for x in pca_transformed], [x[1] for x in pca_transformed], sample_names):
        media_type, genotype, differentiation, additive = metadata[metadata["Name"] == name][["Media", "Cell Line", "Differentiation", "Additives"]].values[0]
        plt.text(x,y,additives_str[additive])
        plt.scatter(x, y, s=40, marker=genotype_marker[genotype], edgecolors=differentation_color[differentiation], facecolor='none' if media_fill[media_type] is None else differentation_color[differentiation])
    explained_vars = model.explained_variance_ratio_
    plt.xlabel("PC 1 " + str(round(explained_vars[0] * 100, 1)) + "%", fontsize=14)
    plt.ylabel("PC 2 " + str(round(explained_vars[1] * 100, 1)) + "%", fontsize=14)
    plt.title(table_name)
    if save_figs:
        plt.savefig("./figures/" + table_name + "_PCA")
    if interactive:
        plt.show()
    plt.clf()

    # do TSNE 
    scaler = StandardScaler()
    transformed = np.log2(feature_table[sample_names].T+1)
    model = TSNE(n_components=2)
    pca_transformed = model.fit_transform(transformed)

    for x,y,name in zip([x[0] for x in pca_transformed], [x[1] for x in pca_transformed], sample_names):
        media_type, genotype, differentiation, additive = metadata[metadata["Name"] == name][["Media", "Cell Line", "Differentiation", "Additives"]].values[0]
        plt.text(x,y,additives_str[additive])
        plt.scatter(x, y, s=40, marker=genotype_marker[genotype], edgecolors=differentation_color[differentiation], facecolor='none' if media_fill[media_type] is None else differentation_color[differentiation])
    plt.title(table_name)
    if save_figs:
        plt.savefig("./figures/" + table_name + "_TSNE")
    if interactive:
        plt.show()
    plt.clf()

significant_features = {}
all_significant_features = []
omnibus_features = []

# do anova
metadata["Genotype"] = metadata["Cell Line"]
metadata["Combined_Media"] = metadata["Media"] + metadata["Additives"]
fields_to_skip = set(feature_table.columns).union(set(metadata.columns))


merged = pd.merge(feature_table.transpose(), metadata, right_on="Name", left_index=True)
print(feature_table.shape)

model_string = '''for_anova ~ C(Genotype) + C(Differentiation) + C(Combined_Media)'''
num_features = feature_table.shape[0]
for column in range(num_features):
    vals = merged.iloc[:,column].values
    new_vals = []
    for x in vals:
        new_vals.append(np.log2(x))
    merged["for_anova"] = new_vals
    model = ols(model_string, data=merged).fit()
    if model.f_pvalue < .05 / (num_features * 7):
        omnibus_features.append(column)
    for index, pvalue in zip(model.pvalues.index[1:], model.pvalues.values[1:]):
        if index not in significant_features:
            significant_features[index] = set()
        if pvalue < .05 / (num_features * 7):
            significant_features[index].add(column)
            all_significant_features.append(column)

differentation_colors = []
genotype_colors = []
media_colors = []
new_names = []
for name in sample_names:
    media_type = metadata[metadata["Name"] == name]["Combined_Media"].values[0]
    if media_type == 'basalNone':
        media_colors.append("forestgreen")
    elif media_type == 'complexNone':
        media_colors.append("firebrick")
    elif media_type == 'basalprolactin':
        media_colors.append("darkorange")
    elif media_type == 'basaldheas':
        media_colors.append("gold")
    genotype = metadata[metadata["Name"] == name]["Cell Line"].values[0]
    if genotype == 'wt':
        genotype_colors.append("aqua")
    elif genotype == 'cyp19a1 KO':
        genotype_colors.append('fuchsia')
    elif genotype == 'prlr KO':
        genotype_colors.append('chartreuse')
    elif genotype == 'ptn KO':
        genotype_colors.append('blue')
    differentiation = metadata[metadata["Name"] == name]["Differentiation"].values[0]
    differentation_colors.append('k' if differentiation == "ipsc" else 'r')
    new_name = ','.join([genotype, differentiation, media_type])
    new_names.append(new_name)

if interactive or save_figs:
    #g = sns.clustermap(new, xticklabels=new_names, yticklabels=master_table.loc[all_significant_features]["annot"], col_colors=[media_colors, genotype_colors, differentation_colors])
    g = sns.clustermap(np.log2(feature_table.iloc[omnibus_features, :][sample_names]), cmap='vlag', xticklabels=new_names, yticklabels=feature_table.iloc[omnibus_features, :]['annot'], z_score=0, col_colors=[media_colors, genotype_colors, differentation_colors])
    g.fig.suptitle(table_name + " all significant OMNIBUS features,\n p-adj(bonferroni) < .05")
    if save_figs:
        plt.savefig("./figures/" + table_name + "_OMNI_clustermap")
    if interactive:
        plt.show()
    plt.clf()

    g = sns.clustermap(np.log2(feature_table.iloc[all_significant_features, :][sample_names]), cmap='vlag', xticklabels=new_names, yticklabels=feature_table.iloc[all_significant_features, :]['annot'], z_score=0, col_colors=[media_colors, genotype_colors, differentation_colors])
    g.fig.suptitle(table_name + " all significant TERM features,\n p-adj(bonferroni) < .05")
    if save_figs:
        plt.savefig("./figures/" + table_name + "_TERM_clustermap")
    if interactive:
        plt.show()
    plt.clf()

groups = [
    # cyp19a1 comparisons
    [("cyp19a1 KO", "ipsc", "basal", "None"), ("wt", "ipsc", "basal", "None")],
    [("cyp19a1 KO", "differentiated", "basal", "None"), ("wt", "differentiated", "basal", "None")],
    [("cyp19a1 KO", "differentiated", "complex", "None"), ("wt", "differentiated", "complex", "None")],
    [("cyp19a1 KO", "differentiated", "basal", "dheas"), ("wt", "differentiated", "basal", "None")],
    [("cyp19a1 KO", "differentiated", "basal", "dheas"), ("cyp19a1 KO", "differentiated", "basal", "None")],

    # prlr comparisons
    [("prlr KO", "ipsc", "basal", "None"), ("wt", "ipsc", "basal", "None")],
    [("prlr KO", "differentiated", "basal", "None"), ("wt", "differentiated", "basal", "None")],
    [("prlr KO", "differentiated", "complex", "None"), ("wt", "differentiated", "complex", "None")],
    [("prlr KO", "differentiated", "basal", "prolactin"), ("wt", "differentiated", "basal", "None")],
    [("prlr KO", "differentiated", "basal", "prolactin"), ("prlr KO", "differentiated", "basal", "None")],

    # ptn comparisons
    [("ptn KO", "ipsc", "basal", "None"), ("wt", "ipsc", "basal", "None")],
    [("ptn KO", "differentiated", "basal", "None"), ("wt", "differentiated", "basal", "None")],
    [("ptn KO", "differentiated", "complex", "None"), ("wt", "differentiated", "complex", "None")],
]    
def ttest(row, indices_a, indices_b):
    #_, p = stats.ttest_ind(row[indices_a], row[indices_b])
    values_A = np.log2([float(x+1) for x in row[indices_a] if not np.isnan(x)])
    values_B = np.log2([float(x+1) for x in row[indices_b] if not np.isnan(x)])
    if len(values_A) + len(values_B) > 3 and list(values_A) and list(values_B):
        t, p = stats.ttest_ind(values_A, values_B)
        if np.isnan(p):
            t, p = 0, 1
        return p
    else:
        return 1

def ttest2(row, indices_a, indices_b):
    #_, p = stats.ttest_ind(row[indices_a], row[indices_b])
    values_A = np.log2([float(x+1) for x in row[indices_a] if not np.isnan(x)])
    values_B = np.log2([float(x+1) for x in row[indices_b] if not np.isnan(x)])
    if len(values_A) + len(values_B) > 3 and list(values_A) and list(values_B):
        t, p = stats.ttest_ind(values_A, values_B)
        if np.isnan(p):
            t, p = 0, 1
        return t
    else:
        return 0
    
def log2fc(row, indices_a, indices_b):
    values_A = np.log2(np.mean(row[indices_a]))
    values_B = np.log2(np.mean(row[indices_b]))
    return values_A - values_B

for group_A_params, group_B_params in groups:
        genotype_A, differentiation_A, media_A, additive_A = group_A_params
        genotype_B, differentiation_B, media_B, additive_B = group_B_params
        #output = open("./analysis/" + ','.join(group_A_params) + "_vs_" + ','.join(group_B_params), 'a+')
        table = feature_table.copy()

        group_A_names = metadata[
                    (metadata["Media"] == media_A) &
                    (metadata["Cell Line"] == genotype_A) &
                    (metadata["Differentiation"] == differentiation_A) &
                    (metadata["Additives"] == additive_A)
                ]["Name"]

        group_B_names = metadata[
                    (metadata["Media"] == media_B) &
                    (metadata["Cell Line"] == genotype_B) &
                    (metadata["Differentiation"] == differentiation_B) &
                    (metadata["Additives"] == additive_B)
                ]["Name"]
        
        group_A = [name for name in group_A_names if name in table.columns]
        group_B = [name for name in group_B_names if name in table.columns]
        names = [name for name in list(group_A_names) + list(group_B_names) if name in table.columns]

        table['p-val'] = table.apply(ttest, axis=1, args=(group_A, group_B))
        table['t-score'] = table.apply(ttest2, axis=1, args=(group_A, group_B))
        table['log2fc'] = table.apply(log2fc, axis=1, args=(group_A, group_B))
        table['significant'], table['corrected'] = multitest.fdrcorrection(table['p-val'], 0.05)

        for_mummichog = pd.DataFrame()
        for_mummichog['m/z'] = table['mz']
        for_mummichog['rtime'] = table['rtime']
        for_mummichog['p-value'] = table['corrected']
        for_mummichog['t-score'] = table["t-score"]
        for_mummichog['custom_id'] = table['annot']
        for_mummichog.to_csv("./mummichog_data/" + table_name + "_".join(group_A_params) + "_vs_" + "_".join(group_B_params) + ".txt", sep="\t", index=False)

        significant_features = pd.DataFrame()
        significant_features['m/z'] = table[table['significant'] == True]['mz']
        significant_features['rtime'] = table[table['significant'] == True]['rtime']
        significant_features['p-value'] = table[table['significant'] == True]['corrected']
        significant_features['t-score'] = table[table['significant'] == True]["t-score"]
        significant_features['custom_id'] = table[table['significant'] == True]['annot']
        significant_features['log2fc'] = table[table['significant'] == True]['log2fc'] 
        significant_features['name'] = table[table['significant'] == True]['annot_names']
        significant_features['emp_ID'] = table[table['significant'] == True]['emp_id']
        if category_class_annots:
            significant_features['categories'] = table[table['significant'] == True]['categories'] 
            significant_features['classes'] = table[table['significant'] == True]['classes'] 
        significant_features.to_csv("./annotated_significant_features/" + table_name + "_".join(group_A_params) + "_vs_" + "_".join(group_B_params) + ".tsv", sep="\t", index=False)

        new_table = pd.DataFrame()
        new_table['m/z'] = table['mz']
        new_table['rtime'] = table['rtime']
        new_table['p-value'] = table['corrected']
        new_table['t-score'] = table["t-score"]
        new_table['custom_id'] = table['annot']
        new_table['log2fc'] = table['log2fc'] 
        new_table['name'] = table['annot_names']
        new_table['emp_ID'] = table['emp_id']
        if category_class_annots:
            new_table['categories'] = table['categories']
            new_table['classes'] = table['classes']
        if category_class_annots:
            fold_change_results = []
            fc_counts = {}
            total_counts = {}
            total_lipids = 0
            total_up = 0
            total_down = 0
            for fc, category_set in zip(new_table['log2fc'], new_table['categories']):
                if type(category_set) is set and len(category_set) == 1:
                    category = list(category_set)[0]
                    if category not in fc_counts:
                        fc_counts[category] = {'up': 0, 'down': 0}
                        total_counts[category] = 0
                    total_counts[category] += 1
                    if fc > 0:
                        total_lipids += 1
                        total_up += 1
                        fc_counts[category]['up'] += 1
                    else:
                        total_lipids += 1
                        total_down += 1
                        fc_counts[category]['down'] += 1
            up_percent = total_up / (total_down + total_up)
            down_percent = 1 - up_percent
            all_pvals = []
            results = []
            for category, up_downs in fc_counts.items():

                observed_up = up_downs['up']
                observed_down = up_downs['down']

                expected_up = total_counts[category] * up_percent
                expected_down = total_counts[category] * down_percent
                fc_counts[category]['expected_up'] = expected_up
                fc_counts[category]['expected_down'] = expected_down
                hpd_down = stats.hypergeom(total_lipids, total_down, observed_up + observed_down)
                hpd_up = stats.hypergeom(total_lipids, total_up, observed_up + observed_down)
                p_val_down = hpd_down.pmf(observed_down)
                p_val_up = hpd_up.pmf(observed_up)
                all_pvals.append(p_val_down)
                all_pvals.append(p_val_up)
                df2 = {'Category': category, 
                       'Expected_Up': expected_up, 
                       'Observed_Up': observed_up, 
                       'Expected_Down': expected_down, 
                       'Observed_Down': observed_down,
                       'pval_up_raw': p_val_up,
                       'pval_down_raw': p_val_down}
                results.append(df2)
            sig, corr_pvals = multitest.fdrcorrection(all_pvals, 0.05)    
            for i, pval in enumerate(corr_pvals):
                r_index = int(i//2)
                if i % 2 == 0:
                    results[r_index]['pval_down_corr'] = pval
                else:
                    results[r_index]['pval_up_corr'] = pval
                    fold_change_results.append(results[r_index])
            fold_change_results = pd.DataFrame(fold_change_results)
            fold_change_results.to_csv("./lipidomics_differential_abundance_analysis/" + table_name + "_".join(group_A_params) + "_vs_" + "_".join(group_B_params) + ".txt", sep="\t", index=False)
