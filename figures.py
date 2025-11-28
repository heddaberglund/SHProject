import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import seaborn as sns

plt.style.use(['science', 'no-latex'])

"""
File for making all plots included in report
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

def plot_necr_rates(necr_rates, lengths):
    plt.plot(lengths, necr_rates)
    plt.title('NECR Rate vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('NECR Rate (Mcps)')
    
    plt.show()
