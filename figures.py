import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import seaborn as sns

plt.style.use(['science', 'no-latex'])

"""
File for making all plots included in report, first six are very standard made using matplotlib, 
plot_seaborn at the end uses seaborn for aesthetics
"""

simulated_counts = 100000  #number of simulated decays for reference line in total counts plot

def plot_site_coords(xis, yis, zis, lengths, actual_site):
    """
    Plot of estimated decay site coordinates (from simple image reconstruction algorithm) vs detector length
    """

    #actual site lines are +-10 in all dir, hashed out lines add shaded area of +-10 round actual site
    #plt.fill_between(lengths, actual_site[0]-10, actual_site[0]+10, color='pink', alpha=0.3)
    #plt.fill_between(lengths, actual_site[1]-10, actual_site[1]+10, color='violet', alpha=0.3)
    #plt.fill_between(lengths, actual_site[2]-10, actual_site[2]+10, color='lightblue', alpha=0.3)
    
    plt.plot(lengths, xis, label='X Coordinate')
    plt.plot(lengths, yis, label='Y Coordinate')
    plt.plot(lengths, zis, label='Z Coordinate')
    plt.axhline(y=actual_site[0], color='pink', linestyle='--', label='Actual X')
    plt.axhline(y=actual_site[1], color='violet', linestyle='--', label='Actual Y')
    plt.axhline(y=actual_site[2], color='lightblue', linestyle='--', label='Actual Z')

    plt.title('Estimated Decay Site Coordinates vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Coordinate Value (cm)')
    plt.legend()
    
    plt.show()

def plot_complex(cs, lengths):
    """
    Complex counts vs. detector length
    """
    plt.plot(lengths, cs)
    plt.title('Events with >2 Hits vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Complex Counts')
    
    plt.show()

def plot_trues(ts, lengths):
    """
    True counts vs. detector length
    """
    plt.plot(lengths, ts)
    plt.title('True Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('True Counts')
    
    plt.show()

def plot_all_counts(ts,ss,rs,cs,lengths):
    """
    Comparison of all our different classifications of counts vs detector length in one plot
    """
    plt.plot(lengths, ts, label='True Counts')
    plt.plot(lengths, ss, label='Scatter Counts')
    plt.plot(lengths, rs, label='Single Counts')
    plt.plot(lengths, cs, label='Complex (>2 hits) Counts')
    plt.title('Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Counts')
    plt.legend()
    
    plt.show()

def plot_total_counts(totals,lengths):
    """
    Total detected counts (irrespective of type) vs detector length
    one loglog plot and one linear plot with simulated counts line for reference
    """

    #loglog plot
    #plt.figure()
    plt.loglog(lengths, totals)
    plt.title('Total Detected Counts vs Detector Length (Log-Log Scale)')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Total Detected Counts')
    plt.grid(True, which="both", ls="--")
    plt.show()
    

    #plinear plot
    plt.axhline(y=simulated_counts, color='r', linestyle='--', label='Simulated Decays')
    plt.plot(lengths, totals, label='Total Detected Counts')
    plt.title('Total Detected Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Total Detected Counts')
    plt.legend()
    
    plt.show()

def plot_necr(ns, lengths):
    
    plt.plot(lengths, ns)
    plt.title('NEC vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('NEC')
    
    plt.show()

def plot_seaborn(lengths, complex, scatter, random, true, necrs, totals, xis, yis, zis, actual_site):
    """
    All the same plots as previously made, but using seaborn for aesthetics in report
    This is admittedly quite repetitive and clunky, but it works for now
    """

    necr_df = pd.DataFrame({"Detector Length (mm)": lengths, "NEC": necrs})
    
    sns.set_theme(style="ticks", context="paper", palette="deep")

    plt.figure(figsize=(6, 4), dpi=300)  

    sns.lineplot(data=necr_df, x="Detector Length (mm)", y="NEC", linewidth=2)

    plt.title("NEC vs Detector Length", fontsize=14, weight="bold")
    plt.xlabel("Detector Length (mm)", fontsize=12)
    plt.ylabel("NEC", fontsize=12)

    sns.despine()             
    plt.tight_layout()        
    plt.show()

    totals_df = pd.DataFrame({"Detector Length (mm)": lengths, "Total Detected Counts": totals})
    plt.figure(figsize=(6, 4), dpi=300)
    #add line at 100000 simulated counts
    plt.axhline(y=simulated_counts, color='r', linestyle='--', label='Simulated Decays')

    sns.lineplot(data=totals_df, x="Detector Length (mm)", y="Total Detected Counts", linewidth=2)

    plt.title("Total Detected Counts vs Detector Length", fontsize=14, weight="bold")
    plt.xlabel("Detector Length (mm)", fontsize=12)
    plt.ylabel("Total Detected Counts", fontsize=12)

    sns.despine()
    plt.tight_layout()
    plt.show()

    trues_df = pd.DataFrame({"Detector Length (mm)": lengths, "True Counts": true})
    plt.figure(figsize=(6, 4), dpi=300)
    sns.lineplot(data=trues_df, x="Detector Length (mm)", y="True Counts", linewidth=2)
    plt.title("True Counts vs Detector Length", fontsize=14, weight="bold")
    plt.xlabel("Detector Length (mm)", fontsize=12)
    plt.ylabel("True Counts", fontsize=12)
    sns.despine()
    plt.tight_layout()
    plt.show()

    allcounts_df = pd.DataFrame({"Detector Length (mm)": lengths,
                                 "True Counts": true,
                                 "Scatter Counts": scatter,
                                 "Random Counts": random,
                                 "Complex Counts": complex})
    melted = pd.melt(
    allcounts_df,
    id_vars=["Detector Length (mm)"],
    var_name="Count Type",
    value_name="Counts"
    )

    plt.figure(figsize=(6, 4), dpi=300)

    sns.lineplot(
        data=melted,
        x="Detector Length (mm)",
        y="Counts",
        hue="Count Type",
        linewidth=2
        )

    plt.title("Counts vs Detector Length", fontsize=14, weight="bold")
    plt.xlabel("Detector Length (mm)", fontsize=12)
    plt.ylabel("Counts", fontsize=12)

    sns.despine()
    plt.tight_layout()
    plt.show()

    #plot decay site
    site_coords_df = pd.DataFrame({"Detector Length (mm)": lengths,
                                   "X Coordinate": xis,
                                   "Y Coordinate": yis,
                                   "Z Coordinate": zis})
    plt.figure(figsize=(6, 4), dpi=300)
    sns.lineplot(data=site_coords_df, x="Detector Length (mm)", y="value", hue="variable", linewidth=2)
    plt.axhline(y=actual_site[0], color='pink', linestyle='--', label='Actual X')
    plt.axhline(y=actual_site[1], color='violet', linestyle='--', label='Actual Y')
    plt.axhline(y=actual_site[2], color='lightblue', linestyle='--', label='Actual Z')
    plt.title("Estimated Decay Site Coordinates vs Detector Length", fontsize=14, weight="bold")
    plt.xlabel("Detector Length (mm)", fontsize=12)
    plt.ylabel("Coordinate Value (cm)", fontsize=12)
    sns.despine()
    plt.tight_layout()
    plt.show()