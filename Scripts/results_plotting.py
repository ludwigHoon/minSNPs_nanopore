# plot result (cc-inference)

from _plotly_utils.colors.cmocean import deep
from _operator import add
from _operator import is_
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import textwrap
# Note: For PNG export, you may need to install: pip install kaleido

def wrap_title(title, width=80):
    """Wrap long titles to multiple lines for better display"""
    return "<br>".join(textwrap.wrap(title, width=width))

short_reads_results = pd.read_csv("../Results/short_reads_individual_results.csv")
simulation_results = pd.read_csv("../Results/simulation_individual_results.csv")
lab_results = pd.read_csv("../Results/lab_individual_results.csv")

results = {
    "short_reads": short_reads_results,
    "simulation": simulation_results,
    "lab": lab_results
}

n_snps_to_set = {
    7: "1 SNP set",
    13: "2 SNP sets",
    18: "3 SNP sets",
    23: "4 SNP sets",
    27: "5 SNP sets",
    99: "30 SNP sets",
    199: "79 SNP sets",
    300: "147 SNP sets",
    385: "200 SNP sets"
}

# Harmonise column names for plotting
# random ->  random_400
# minsnps -> non_random_200
short_reads_results.loc[:, "snps_type"] = short_reads_results.loc[:, "search_type"].map({
    "random": "random_400",
    "minsnps": "non_random_200"})

def transform_snps_number(row):
    if row["snps_type"] == "random_400":
        return row["n_snps"]
    elif row["snps_type"] == "non_random_200":
        mapping = {
            1: 7,
            2: 13,
            3: 18,
            4: 23,
            5: 27,
            30: 99,
            79: 199,
            147: 300,
            200: 385
        }
        return mapping.get(row["n_snps"])
    else:
        raise ValueError(f"Unknown snps_type: {row['snps_type']}")

short_reads_results.loc[:, "snps_used"] = short_reads_results.apply(transform_snps_number, axis=1)
simulation_results.loc[:, "sample"] = simulation_results.loc[:, "pubmlst_id"]
lab_results.loc[:, "sample"] = lab_results.loc[:, "barcode_id"]
lab_results.loc[:, "CC_filled"] = lab_results.loc[:, "CC"]
lab_results.loc[:, "n_reads"] = lab_results.loc[:, "max_reads_used"]
simulation_results.loc[:, "n_reads"] = simulation_results.loc[:, "reads_used"]
lab_results.loc[:, "basecall_method"] = lab_results.loc[:, "call_type"]

def wrap_title(title, width=80):
    """Wrap long titles to multiple lines for better display"""
    return "<br>".join(textwrap.wrap(title, width=width))

def plot_correctness_result(results, title, additional_filtering, output_file):
    df = results.copy()
    if additional_filtering:
        df = df.query(additional_filtering)
    
    # Correctness = number of samples where CC_filled == predicted_MLST.CC
    df['is_correct'] = df['CC_filled'] == df['predicted_MLST.CC']
    
    df_summarised = df.groupby(["n_reads", "snps_type", "snps_used"]).agg({
        'is_correct': 'sum',
        'sample': 'nunique',  # Count unique samples instead of listing them
    }).reset_index()

    df_summarised['proportion'] = (df_summarised['is_correct'] / df_summarised['sample']) * 100
    # Flatten column names
    df_summarised.columns = ['n_reads', 'snps_type', 'snps_used', 'correct_predictions', 'total_samples', 'proportion']
    df_summarised["n_snps_cat"] = df_summarised["snps_used"].apply(lambda x: 100 if x == 99 else 200 if x == 199 else 400 if x == 385 else x)
    df_summarised['snp_label'] = df_summarised.apply(
        lambda x: f"{x['snps_used']} SNPs ({n_snps_to_set.get(x['snps_used'], '')})" if x['snps_type'] == 'non_random_200' else f"{x['n_snps_cat']} SNPs",
        axis=1
    )
    # Create color mapping for different numbers of SNPs
    unique_snps = sorted(df_summarised['n_snps_cat'].unique())
    color_palette = px.colors.qualitative.Plotly
    snp_color_map = {snp: color_palette[i % len(color_palette)] for i, snp in enumerate(unique_snps)}

    # Calculate out-performance: non_random_200 - random_400
    df_outperformance = df_summarised.copy(deep=True)
    
    # Pivot to get both SNP types for comparison
    df_pivot = df_outperformance.pivot_table(
        index=['n_reads', 'n_snps_cat'], 
        columns='snps_type', 
        values='correct_predictions', 
        fill_value=0
    ).reset_index()
    
    # Calculate outperformance (non_random_200 - random_400)
    if 'non_random_200' in df_pivot.columns and 'random_400' in df_pivot.columns:
        df_pivot['outperformance'] = df_pivot['non_random_200'] - df_pivot['random_400']
    else:
        df_pivot['outperformance'] = 0
    
    fig = go.Figure()
    
    for n_snps_cat in df_pivot['n_snps_cat'].unique():
        data_subset = df_pivot[df_pivot['n_snps_cat'] == n_snps_cat]
        snp_color = snp_color_map[n_snps_cat]
        
        fig.add_trace(
            go.Bar(
                x=data_subset['n_reads'],
                y=data_subset['outperformance'],
                name=f'{"99/100" if n_snps_cat == 100 else "199/200" if n_snps_cat == 200 else "385/400" if n_snps_cat == 400 else n_snps_cat} SNPs',
                marker_color=snp_color,
            )
        )

    wrapped_title = wrap_title(f"{title.replace('Number of samples assigned to correct lineage', 'Outperformance of resolution-optimised against random SNPs')}", width=70)

    fig.update_layout(
        title=wrapped_title,
        xaxis_title="Number of reads used",
        yaxis_title="Outperformance (Number of samples)",
        legend_title_text="Number of SNPs used",
        width=800,
        height=600,
        font=dict(size=12),
        margin=dict(b=80, r=150)  # Add bottom margin for rotated labels, right margin for legend
    )
    
    # Ensure x-axis labels are visible and not shrunk
    fig.update_xaxes(
        tickmode='linear',
        tick0=1000,
        dtick=1000,  # Show every 1000 reads
        title_standoff=25  # Add space between axis and title
    )
    
    # Position legend
    fig.update_layout(legend_x=1.05, legend_y=1.0)
    fig.update_layout(legend_xanchor="left", legend_yanchor="top")
    
    # Save files
    outperformance_html = output_file.replace('.html', '_outperformance.html')
    outperformance_png = output_file.replace('.html', '_outperformance.png')
    
    fig.write_html(outperformance_html)
    fig.write_image(
        outperformance_png,
        width=800,
        height=600,
        scale=3,
        format='png'
    )
    
    print(f"Outperformance plot saved to {outperformance_html}")
    print(f"Outperformance PNG saved to {outperformance_png}")


    for snp_type in df_summarised['snps_type'].unique():
    # Create subplot with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        # Add bar chart for correct predictions (left y-axis)
        for n_snps in df_summarised['snps_used'].unique():
            data_subset = df_summarised[
                (df_summarised['snps_type'] == snp_type) &
                (df_summarised['snps_used'] == n_snps)
            ]

            if data_subset.empty:
                continue
            # Use mapped color for this specific number of SNPs
            snp_color = snp_color_map[data_subset['n_snps_cat'].values[0]]

            fig.add_trace(
                go.Bar(
                    x=data_subset['n_reads'],
                    y=data_subset['correct_predictions'],
                    name=data_subset['snp_label'].values[0],
                    legendgroup=f"{n_snps}",
                    marker_color=snp_color,
                    opacity=0.7
                ),
                secondary_y=False,
            )
            
            # Add line chart for proportion (right y-axis) - no legend since same color as bar
            # fig.add_trace(
            #     go.Scatter(
            #         x=data_subset['n_reads'],
            #         y=data_subset['proportion'],
            #         mode='lines+markers',
            #         name=f'{n_snps} SNPs - %',
            #         legendgroup=f"{n_snps}",
            #         showlegend=False,  # Hide from legend since bar already shows this
            #         line_color=snp_color,
            #         line_width=3,
            #         marker_color=snp_color
            #     ),
            #     secondary_y=True,
            # )
        
        # Set x-axis title
        fig.update_xaxes(
            title_text="Number of reads used",
            tickmode='linear',
            tick0=1000,
            dtick=1000,  # Show every 1000 reads
            title_standoff=25  # Add space between axis and title
        )
        
        # add h-line 
        fig.add_hline(
            y=df_summarised["total_samples"].max(),  # Horizontal line at 50%
            line_dash="dash",
            line_color="red",
            annotation_text="Number of tested samples",
            annotation_position="top left"
        )

        # Set y-axes titles
        fig.update_yaxes(title_text="Number of correctly assigned samples", secondary_y=False)
        #fig.update_yaxes(title_text="Percentage of correctly assigned samples (%)", secondary_y=True)
        
        # Update layout
        snp_cat = {
            "random_400": "random SNPs",
            "non_random_200": "resolution-optimised SNPs"
        }
        full_title = f"{title} using {snp_cat[snp_type]}"
        wrapped_title = wrap_title(full_title, width=70)  # Adjust width for 800px figure
        
        fig.update_layout(
            title=wrapped_title,
            xaxis_title="Number of reads used",
            legend_title_text="Number of SNPs used",
            width=800,  # Set figure width for publication
            height=600,  # Set figure height for publication
            font=dict(size=12),  # Set font size for readability
            margin=dict(b=80, r=150)  # Add bottom margin for rotated labels, right margin for legend
        )
        
        # Position legend with gap from secondary y-axis
        fig.update_layout(legend_x=1.05, legend_y=1.0)
        fig.update_layout(legend_xanchor="left", legend_yanchor="top")
        
        # Save as HTML
        fig.write_html(f"{snp_type}_{output_file}")
        
        # Save as PNG with publication quality settings
        png_filename = f"{snp_type}_{output_file}".replace('.html', '.png')
        fig.write_image(
            png_filename,
            width=800,
            height=600,
            scale=3,  # Scale factor for higher resolution (3x = 300 DPI equivalent)
            format='png'
        )
        
        print(f"Plot saved to {snp_type}_{output_file}")
        print(f"PNG saved to {png_filename}")


for experiment in results.keys():
    if experiment == "lab":
        additional_filtering = "call_type == 'fast'"
    else:
        additional_filtering = None
    titles = {
        "short_reads": "Short-reads",
        "simulation": "Simulated long reads",
        "lab": "Lab Results (fast basecalling)"
    }
    experiment_title = titles[experiment]
    
    # Generate correctness plots
    plot_correctness_result(
        results[experiment],
        title=f"{experiment_title} - Number of samples assigned to correct lineage vs number of reads",
        additional_filtering=additional_filtering,
        output_file=f"{experiment}_correctness_plot.html"
    )



def plot_correctness_result_combined(results, title, additional_filtering=None):
    combined_fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=[
            "A. Fast basecalling",
            "B. High accuracy (hac) basecalling", 
            "C. Super high accuracy (sup) basecalling"
        ],
        shared_yaxes=True,
        vertical_spacing=0.08
    )
    df = results.copy()
    if additional_filtering:
        df = df.query(additional_filtering)
    
    # Correctness = number of samples where CC_filled == predicted_MLST.CC
    df['is_correct'] = df['CC_filled'] == df['predicted_MLST.CC']
    
    df_summarised = df.groupby(["basecall_method", "n_reads", "snps_type", "snps_used"]).agg({
        'is_correct': 'sum',
        'sample': 'nunique',  # Count unique samples instead of listing them
    }).reset_index()

    df_summarised['proportion'] = (df_summarised['is_correct'] / df_summarised['sample']) * 100
    # Flatten column names
    df_summarised.columns = ["basecall_method", 'n_reads', 'snps_type', 'snps_used', 'correct_predictions', 'total_samples', 'proportion']
    df_summarised["n_snps_cat"] = df_summarised["snps_used"].apply(lambda x: 100 if x == 99 else 200 if x == 199 else 400 if x == 385 else x)
    df_summarised['label'] = df_summarised.apply(
        lambda x: f"{x['snps_used']} SNPs ({n_snps_to_set.get(x['snps_used'], '')})" if x['snps_type'] == 'non_random_200' else f"{x['n_snps_cat']} SNPs",
        axis=1
    )
    # Create color mapping for different numbers of SNPs
    unique_snps = sorted(df_summarised['n_snps_cat'].unique())
    color_palette = px.colors.qualitative.Plotly
    snp_color_map = {snp: color_palette[i % len(color_palette)] for i, snp in enumerate(unique_snps)}

    # Map basecall methods to subplot rows
    basecall_row_map = {
        'fast': 1,
        'hac': 2,
        'sup': 3
    }


    for basecall_method in df_summarised['basecall_method'].unique():
        if basecall_method not in basecall_row_map:
            continue
            
        row_num = basecall_row_map[basecall_method]
        
        # Filter data for this SNP type and basecall method
        data_subset = df_summarised[
            (df_summarised['basecall_method'] == basecall_method)
        ]
        
        if data_subset.empty:
            continue
            
        # Group by SNP count for this combination
        for n_snps_cat in data_subset['n_snps_cat'].unique():
            snp_subset = data_subset[data_subset['n_snps_cat'] == n_snps_cat]
            snp_color = snp_color_map[n_snps_cat]
            
            # Create legend label with SNP type info
            legend_label = snp_subset['label'].values[0]

            combined_fig.add_trace(
                go.Bar(
                    x=snp_subset['n_reads'],
                    y=snp_subset['correct_predictions'],
                    name=legend_label,
                    marker_color=snp_color,
                    opacity=0.7,
                    showlegend=(row_num == 2),  # Only show legend for middle subplot
                    legendgroup=f"{n_snps_cat}"
                ),
                row=row_num, col=1
            )

    # Update layout
    combined_fig.update_layout(
        title=wrap_title(title, width=90),
        width=900,
        height=1000,
        font=dict(size=11),
        margin=dict(b=80, r=150, l=80),
        showlegend=True
    )

    # Update x-axes
    combined_fig.update_xaxes(
        title_text="Number of reads used",
        tickmode='linear',
        tick0=0,
        dtick=1000,
        tickangle=0,
        title_standoff=25
    )

    # Update y-axes
    combined_fig.update_yaxes(
        title_text="Number of correctly assigned samples",
        row=2, col=1  # Only set y-axis title for middle subplot
    )

    # Add horizontal line for total samples on each subplot
    for row in range(1, 4):
        max_samples = df_summarised["total_samples"].max()
        combined_fig.add_hline(
            y=max_samples,
            line_dash="dash",
            line_color="red",
            annotation_text="Total samples tested",
            annotation_position="top left",
            row=row, col=1
        )

    # Position legend to the right and middle
    combined_fig.update_layout(
        legend=dict(
            title_text="Number of SNPs used",
            x=1.02,
            y=0.5,
            xanchor="left",
            yanchor="middle"
        )
    )

    # Save combined figure
    output_html = "lab_correctness_combined.html"
    output_png = "lab_correctness_combined.png"
    
    combined_fig.write_html(output_html)
    combined_fig.write_image(
        output_png,
        width=900,
        height=1200,
        scale=3,
        format='png'
    )

    print(f"Combined correctness plot saved to {output_html}")
    print(f"Combined PNG saved to {output_png}")


# Generate the combined correctness plot for lab results
plot_correctness_result_combined(
    results['lab'],
    title="Lab-generated long reads - Number of samples assigned to correct lineage vs number of reads using resolution-optimised SNPs",
    additional_filtering="snps_type == 'non_random_200'"
)






#  Gene results
gene_short_reads_results = pd.read_csv("../Results/short_reads_gene_results_individual.csv")
gene_simulation_results = pd.read_csv("../Results/simulated_gene_result.csv")
gene_lab_results = pd.read_csv("../Results/lab_gene_result.csv")

gene_results = {
    "short_reads": gene_short_reads_results,
    "simulation": gene_simulation_results,
    "lab": gene_lab_results
}

# Harmonise column names for plotting
gene_simulation_results.loc[:, "sample"] = gene_simulation_results.loc[:, "seq_name"]
gene_lab_results.loc[:, "sample"] = gene_lab_results.loc[:, "Barcode.ID"]
gene_simulation_results.loc[:, "n_reads"] = gene_simulation_results.loc[:, "reads_used"]
gene_lab_results.loc[:, "n_reads"] = gene_lab_results.loc[:, "max_reads_used"]
gene_simulation_results.loc[:, "n_kmer"] = gene_simulation_results.loc[:, "n_sequence_used"]
gene_lab_results.loc[:, "n_kmer"] = gene_lab_results.loc[:, "n_sequence_used"]

# ACTUAL
# mecA(AMRFINDER),lukF.PV(AMRFINDER),lukS.PV(AMRFINDER)
gene_lab_results.loc[:, "mecA(AMRFINDER)"] = gene_lab_results.loc[:, "MRSA"]
gene_lab_results.loc[:, "lukF.PV(AMRFINDER)"] = gene_lab_results.loc[:, "PVL"]
gene_lab_results.loc[:, "lukS.PV(AMRFINDER)"] = gene_lab_results.loc[:, "PVL"]
gene_simulation_results.loc[:, "mecA(AMRFINDER)"] = gene_simulation_results.loc[:, "mecA"]
gene_simulation_results.loc[:, "lukF.PV(AMRFINDER)"] = gene_simulation_results.loc[:, "lukF_PV"]
gene_simulation_results.loc[:, "lukS.PV(AMRFINDER)"] = gene_simulation_results.loc[:, "lukS_PV"]

# FOUND
# lukF_PV,lukS_PV,mecA

gene_simulation_results.loc[:, "lukF_PV"] = gene_simulation_results.loc[:, "lukF_prop"]
gene_simulation_results.loc[:, "lukS_PV"] = gene_simulation_results.loc[:, "lukS_prop"]
gene_simulation_results.loc[:, "mecA"] = gene_simulation_results.loc[:, "mecA_prop"]

gene_lab_results.loc[:, "lukF_PV"] = gene_lab_results.loc[:, "lukF_prop"]
gene_lab_results.loc[:, "lukS_PV"] = gene_lab_results.loc[:, "lukS_prop"]
gene_lab_results.loc[:, "mecA"] = gene_lab_results.loc[:, "mecA_prop"]



def plot_gene_sensitivity_specificity(
    results, 
    title, 
    additional_filtering=None,
    detection_threshold=0.6,
    output_file="gene_sensitivity_plot.html"
):
    df = results.copy()
    if additional_filtering:
        df = df.query(additional_filtering)
    
    genes_predicted_actual = {
        'mecA': 'mecA(AMRFINDER)',
        'lukF_PV': 'lukF.PV(AMRFINDER)',
        'lukS_PV': 'lukS.PV(AMRFINDER)'
    }

    # Create binary columns for each gene based on detection threshold
    for gene in genes_predicted_actual.keys():
        df[f'{gene}_detected'] = df[gene] >= detection_threshold

    # Calculate TP, TN, FP, FN for each gene
    for gene_pred, gene_actual in genes_predicted_actual.items():
        # True Positives: predicted positive AND actual positive
        df[f'{gene_pred}_TP'] = (df[f'{gene_pred}_detected'] == True) & (df[gene_actual] == True)
        
        # True Negatives: predicted negative AND actual negative  
        df[f'{gene_pred}_TN'] = (df[f'{gene_pred}_detected'] == False) & (df[gene_actual] == False)
        
        # False Positives: predicted positive BUT actual negative
        df[f'{gene_pred}_FP'] = (df[f'{gene_pred}_detected'] == True) & (df[gene_actual] == False)
        
        # False Negatives: predicted negative BUT actual positive
        df[f'{gene_pred}_FN'] = (df[f'{gene_pred}_detected'] == False) & (df[gene_actual] == True)

    agg_dict = {'sample': 'nunique'}  # Count unique samples    
    # Add TP, TN, FP, FN for each gene to aggregation
    for gene_pred in genes_predicted_actual.keys():
        agg_dict[f'{gene_pred}_TP'] = 'sum'
        agg_dict[f'{gene_pred}_TN'] = 'sum' 
        agg_dict[f'{gene_pred}_FP'] = 'sum'
        agg_dict[f'{gene_pred}_FN'] = 'sum'
    
    df_summary = df.groupby(['n_kmer', 'n_reads']).agg(agg_dict).reset_index()

    # Calculate sensitivity and specificity for each gene in df_summary
    for gene_pred in genes_predicted_actual.keys():
        # Sensitivity = TP / (TP + FN)
        df_summary[f'{gene_pred}_sensitivity'] = (
            df_summary[f'{gene_pred}_TP'] / 
            (df_summary[f'{gene_pred}_TP'] + df_summary[f'{gene_pred}_FN'])
        ).fillna(0)  # Handle division by zero
        
        # Specificity = TN / (TN + FP)
        df_summary[f'{gene_pred}_specificity'] = (
            df_summary[f'{gene_pred}_TN'] / 
            (df_summary[f'{gene_pred}_TN'] + df_summary[f'{gene_pred}_FP'])
        ).fillna(0)  # Handle division by zero
    
    # Define colors for each gene (darker for sensitivity, lighter for specificity)
    gene_colors = {
        'mecA': {'dark': 'darkblue', 'light': 'lightblue'},
        'lukF_PV': {'dark': 'darkgreen', 'light': 'lightgreen'},
        'lukS_PV': {'dark': 'darkred', 'light': 'lightcoral'}
    }
    
    # Create a single combined plot with all genes
    fig = go.Figure()
    
    for gene_pred in genes_predicted_actual.keys():
        colors = gene_colors[gene_pred]
        
        # Add sensitivity line (darker color, solid line)
        fig.add_trace(go.Scatter(
            x=df_summary['n_reads'],
            y=df_summary[f'{gene_pred}_sensitivity'],
            mode='lines+markers',
            name=f'{gene_pred} Sensitivity',
            line_color=colors['dark'],
            line_width=3,
            marker_size=8,
            line_dash='solid'
        ))
        
        # Add specificity line (lighter color, dotted line)
        fig.add_trace(go.Scatter(
            x=df_summary['n_reads'],
            y=df_summary[f'{gene_pred}_specificity'],
            mode='lines+markers',
            name=f'{gene_pred} Specificity',
            line_color=colors['light'],
            line_width=3,
            marker_size=8,
            line_dash='dot'
        ))
    wrapped_title = wrap_title(f"{title.replace(')', '')}, threshold = {detection_threshold:.2f})", width=60)

    fig.update_layout(
        title=wrapped_title,
        xaxis_title="Number of reads used",
        yaxis_title="Sensitivity / Specificity",
        legend_title_text="",
        width=800,
        height=600,
        font=dict(size=12),
        margin=dict(b=80, r=150),
        yaxis=dict(range=[0, 1])  # Set y-axis range from 0 to 1
    )
    
    # Update x-axis formatting
    fig.update_xaxes(
        tickmode='linear',
        tick0=0,
        dtick=1000,
        tickangle=45,
        title_standoff=25
    )
    
    # Position legend
    fig.update_layout(legend_x=1.05, legend_y=1.0)
    fig.update_layout(legend_xanchor="left", legend_yanchor="top")
    
    # Save files
    fig.write_html(output_file)
    fig.write_image(
        output_file.replace('.html', '.png'),
        width=800,
        height=600,
        scale=3,
        format='png'
    )
    
    print(f"Gene performance plot saved to {output_file}")
    print(f"PNG saved to {output_file.replace('.html', '.png')}")
    
    return df_summary, fig

data_sim_k100_t0_6, fig_sim_k100_t0_6 = plot_gene_sensitivity_specificity(
    gene_results['simulation'],
    title="Simulated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = 100)",
    additional_filtering="n_kmer == '100'",
    detection_threshold=0.6,
    output_file="sim_gene_nkmer_100_thres_0.6.html"
)

data_sim_kall_t0_6, fig_sim_kall_t0_6 = plot_gene_sensitivity_specificity(
    gene_results['simulation'],
    title="Simulated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = all)",
    additional_filtering="n_kmer == 'all'",
    detection_threshold=0.6,
    output_file="sim_gene_nkmer_all_thres_0.6.html"
)

data_sim_k100_t0_7, fig_sim_k100_t0_7 = plot_gene_sensitivity_specificity(
    gene_results['simulation'],
    title="Simulated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = 100)",
    additional_filtering="n_kmer == '100'",
    detection_threshold=0.7,
    output_file="sim_gene_nkmer_100_thres_0.7.html"
)

data_lab_k100_t0_6, fig_lab_k100_t0_6 = plot_gene_sensitivity_specificity(
    gene_results['lab'],
    title="Lab-generated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = 100)",
    additional_filtering="n_kmer == '100' & basecall_method == 'fast'",
    detection_threshold=0.6,
    output_file="lab_gene_nkmer_100_thres_0.6.html"
)

data_lab_kall_t0_6, fig_lab_kall_t0_6 = plot_gene_sensitivity_specificity(
    gene_results['lab'],
    title="Lab-generated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = all)",
    additional_filtering="n_kmer == 'all' & basecall_method == 'fast'",
    detection_threshold=0.6,
    output_file="lab_gene_nkmer_all_thres_0.6.html"
)

data_lab_k100_t0_7, fig_lab_k100_t0_7 = plot_gene_sensitivity_specificity(
    gene_results['lab'],
    title="Lab-generated long-read - Gene detection sensitivity/specificity vs number of reads (kmer-used = 100)",
    additional_filtering="n_kmer == '100' & basecall_method == 'fast'",
    detection_threshold=0.7,
    output_file="lab_gene_nkmer_100_thres_0.7.html"
)


# Create combined subplot figure
from plotly.subplots import make_subplots

# Define colors for each gene (consistent across all subplots)
gene_colors = {
    'mecA': {'dark': 'darkblue', 'light': 'lightblue'},
    'lukF_PV': {'dark': 'darkgreen', 'light': 'lightgreen'},
    'lukS_PV': {'dark': 'darkred', 'light': 'lightcoral'}
}

genes_predicted_actual = {
    'mecA': 'mecA(AMRFINDER)',
    'lukF_PV': 'lukF.PV(AMRFINDER)',
    'lukS_PV': 'lukS.PV(AMRFINDER)'
}


datasets_sim = [
    (data_sim_k100_t0_6, 1),
    (data_sim_kall_t0_6, 2), 
    (data_sim_k100_t0_7, 3)
]

datasets_lab = [
    (data_lab_k100_t0_6, 1),
    (data_lab_kall_t0_6, 2),
    (data_lab_k100_t0_7, 3)
]


thresholds = [a/100 for a in range(5, 61, 5)]
datasets_sr = []
for threshold in thresholds:
    data, _ = plot_gene_sensitivity_specificity(
        gene_results['short_reads'],
        title="Tanzania short-reads - Gene detection sensitivity/specificity vs number of reads (kmer-used = 100)",
        additional_filtering="n_kmer == '100'",
        detection_threshold=threshold,
        output_file=f"sr_gene_nkmer_100_thres_{threshold}.html"
    )
    datasets_sr.append((threshold, data))

def combined_gene_plot(datasets, genes_predicted_actual, gene_colors, title, output_file, subplot_titles=None):
    # Create subplot figure with 3 rows and 1 column
    combined_fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=subplot_titles or [
            "A. k-mer = 100, threshold = 0.6",
            "B. k-mer = all, threshold = 0.6", 
            "C. k-mer = 100, threshold = 0.7"
        ],
        shared_yaxes=True,
        horizontal_spacing=0.08
    )

    for data, row_num in datasets:
        for gene_pred in genes_predicted_actual.keys():
            colors = gene_colors[gene_pred]
            
            # Add sensitivity line (darker color, solid line)
            combined_fig.add_trace(
                go.Scatter(
                    x=data['n_reads'],
                    y=data[f'{gene_pred}_sensitivity'],
                    mode='lines+markers',
                    name=f'{gene_pred} Sensitivity',
                    line_color=colors['dark'],
                    line_width=2,
                    marker_size=6,
                    line_dash='solid',
                    showlegend=(row_num == 2),  # Only show legend for first subplot
                    legendgroup=f'{gene_pred}_sens'
                ),
                row=row_num, col=1
            )
            
            # Add specificity line (lighter color, dotted line)
            combined_fig.add_trace(
                go.Scatter(
                    x=data['n_reads'],
                    y=data[f'{gene_pred}_specificity'],
                    mode='lines+markers',
                    name=f'{gene_pred} Specificity',
                    line_color=colors['light'],
                    line_width=2,
                    marker_size=6,
                    line_dash='dot',
                    showlegend=(row_num == 2),  # Only show legend for first subplot
                    legendgroup=f'{gene_pred}_spec'
                ),
                row=row_num, col=1
            )

    # Update layout
    combined_fig.update_layout(
        title=wrap_title(f"{title} - Gene detection: sensitivity/specificity vs number of reads", width=80),
        width=1400,  # Wider to accommodate 3 subplots
        height=500,
        font=dict(size=11),
        margin=dict(b=80, r=150, l=80),
        showlegend=True
    )

    # Update x-axes
    combined_fig.update_xaxes(
        title_text="Number of reads used",
        tickmode='linear',
        tick0=0,
        dtick=1000,  # Larger intervals for smaller subplots
        tickangle=45,
        title_standoff=25
    )

    # Update y-axes
    combined_fig.update_yaxes(
        title_text="Sensitivity / Specificity",
        range=[0, 1],
        row=2, col=1 
    )

    # Update all y-axes ranges
    for row in range(1, 4):
        combined_fig.update_yaxes(range=[0, 1], row=row, col=1)

    # Position legend
    combined_fig.update_layout(
        legend=dict(
            x=1.02,
            y=0.5,
            xanchor="left",
            yanchor="middle"
        )
    )

    # Save combined figure
    combined_fig.write_html(output_file)
    combined_fig.write_image(
        output_file.replace('.html', '.png'),
        width=800,
        height=1200,
        scale=3,
        format='png'
    )

    print(f"Combined gene performance plot saved to {output_file}")
    print(f"Combined PNG saved to {output_file.replace('.html', '.png')}")

combined_gene_plot(datasets_sim, genes_predicted_actual, gene_colors, "Simulated long-read", "sim_gene_combined_performance.html")
combined_gene_plot(datasets_lab, genes_predicted_actual, gene_colors, "Lab-generated long-read", "lab_gene_combined_performance.html")
subplot_titles = [
    "A. k-mer = 100, threshold = 0.05",
    "B. k-mer = 100, threshold = 0.1",
    "C. k-mer = 100, threshold = 0.15"
]
to_plot = [(dat, i) for i, dat in enumerate([data for threshold, data in datasets_sr if threshold in [0.05, 0.1, 0.15]])]
combined_gene_plot(to_plot, genes_predicted_actual, gene_colors, "Short-read", "sr_gene_combined_performance.html", subplot_titles=subplot_titles)