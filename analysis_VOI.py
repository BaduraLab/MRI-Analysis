# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
reference_structure_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail.csv')
structure = pd.read_csv(reference_structure_path)


# Plot boxplot plots for VOIs


volume_name = 'Lobule II'
ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
plt.ylabel('$mm^3$')
plt.xlabel('Genotype')
plt.title(volume_name + ' volumes')
plt.suptitle('') # that's what you're after
# ax.set_xticklabels(['WT', 'KO'])
# plt.show()
plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))

volume_name = 'Substantia nigra, compact part'
ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
plt.ylabel('$mm^3$')
plt.xlabel('Genotype')
plt.title(volume_name + ' volumes')
plt.suptitle('') # that's what you're after
# ax.set_xticklabels(['WT', 'KO'])
# plt.show()
plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))

volume_name = 'Substantia nigra, reticular part'
ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
plt.ylabel('$mm^3$')
plt.xlabel('Genotype')
plt.title(volume_name + ' volumes')
plt.suptitle('') # that's what you're after
# ax.set_xticklabels(['WT', 'KO'])
# plt.show()
plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))

# Plotting by genotype and sex
volume_name = 'Lobules IV-V'
ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
plt.ylabel('$mm^3$')
plt.xlabel('Genotype and Sex')
plt.title(volume_name + ' volumes')
plt.suptitle('') # that's what you're after
# ax.set_xticklabels(['WT', 'KO'])
# plt.show()
plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotypeSex'))