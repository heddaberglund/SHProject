import numpy as np
import pandas as pd

from figures import *

EVENT_ID, GEOMETRY_ID, ENERGY, TIME, R, PHI, Z = 0, 1, 2, 3, 4, 5, 6 #units: _, _, keV, ns, mm, rad
detector_lengths_mm = [100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050, 1100,1150,1200,1250,1300,1350,1400]
files = ["hits100mm", "hits150mm", "hits200mm", "hits250mm", "hits300mm", \
"hits350mm", "hits400mm", "hits450mm", "hits500mm", "hits550mm", "hits600mm", "hits650mm", "hits700mm", "hits750mm", \
"hits800mm", "hits850mm", "hits900mm", "hits950mm", "hits1000mm", "hits1050mm", "hits1100mm", "hits1150mm", "hits1200mm", \
"hits1250mm", "hits1300mm", "hits1350mm", "hits1400mm"]
actual_site = (200.3/4, 100.0/2, 0.0)

def find_true(event):
    """Determine if an event is a true coincidence.
    takes touple of two detections (event) and checks if both energies are within 70kev of 511kev.
    possibility of adding angle check as well, but left out because we don't want to preemptively filter
    out events we use for image reconstruction later on. also simulated tumour signal is not centered on the origin
    in sim so adding a 180 degree angle filter would be inaccurate as true counts do not necessarily produce 180degree angle.
    """
    energy1 = float(event[0][ENERGY])
    energy2 = float(event[1][ENERGY])
    angle = find_angle(event)
    
    #energy filter relatively lenient at 70keV, could be made smaller
    is_true = (abs(energy1 - 511) < 70) and (abs(energy2 - 511) < 70) #and (abs(abs(angle) - 180) < 50)
    return is_true

def find_angle(event):
    """
    simple trig to find angle between two hits in an event based on r, phi coords 
    not used in current analysis but could be useful for future work
    """

    phi0 = float(event[0][PHI])
    phi1 = float(event[1][PHI])
    r0 = float(event[0][R])
    r1 = float(event[1][R])

    x1 = r0*np.cos(phi0)
    y1 = r0*np.sin(phi0)
    
    x2 = r1*np.cos(phi1)
    y2 = r1*np.sin(phi1)
    
    delta_x = x2 - x1
    delta_y = y2 - y1
    
    angle = np.arctan2(delta_y, delta_x) * (180 / np.pi)  # convert to degrees
    return angle

def find_necr(t,s,r,k):
    """Calculate NECR given true (t), scatter (s), random (r) counts, and scale factor k."""
    return (t*t) / (t + s + k*r)

def necr_as_rate(necr, dose, simulated_counts):
    #rate in MBq
    rate = necr * dose / simulated_counts
    return rate

def find_LOR(event):
    """
    finds line of response from two hits, very similar to find_angle so could be made into 
    one fxn returning angle and coords...
    """
    phi0 = float(event[0][PHI])
    phi1 = float(event[1][PHI])
    r0 = float(event[0][R])
    r1 = float(event[1][R])

    x1 = r0*np.cos(phi0)
    y1 = r0*np.sin(phi0)
    
    x2 = r1*np.cos(phi1)
    y2 = r1*np.sin(phi1)

    z1 = float(event[0][Z])
    z2 = float(event[1][Z])

    return (x1, y1, z1), (x2, y2, z2)

def find_midpoint(line):
    """
    finds midpoint of an LOR
    """
    (x1, y1, z1), (x2, y2, z2) = line
    x_mid = (x1 + x2) / 2
    y_mid = (y1 + y2) / 2
    z_mid = (z1 + z2) / 2
    return (x_mid, y_mid, z_mid)

def find_decay_site(midpoints):
    """
    averages midpoints to estimate decay site.. very much a simplification
    """
    n = len(midpoints)
    x_sum = 0
    y_sum = 0
    z_sum = 0
    for midpoint in midpoints:
        x_sum += midpoint[0]
        y_sum += midpoint[1]
        z_sum += midpoint[2]
    return (x_sum / n, y_sum / n, z_sum / n)

def compute_decay_site_difference(estimated, actual):
    "compute difference between estimated and actual decay site"
    dx = abs(estimated[0] - actual[0])
    dy = abs(estimated[1] - actual[1])
    dz = abs(estimated[2] - actual[2])
    return (dx, dy, dz)

def group_counts(df):
    """This fxn groups counts by event ID, and determine whether each event is
    - true: if event id has 2 hits with energy close to 511keV, or event id 
            has 1 hit and next event id has 1 hit and both hits have energy close to 511keV
    - scatter: event id has one hit, or two hits but one or both energies outside 511keV
    - random: event id has 1 hit inconsistent with 511keV or next event id has >1 hit
    - complex: event id has >2 hits... in future could try to classify some of these as true/scatter as complex
               are currently not used for necr stats
    """
    
    LOR_list = []
    event_dict = {}
    for detection in df:
        hit = detection[0].split()
        event_id = hit[EVENT_ID]
        if event_id not in event_dict:
            event_dict[event_id] = []
        event_dict[event_id].append(hit)
    
    complex_no = 0
    group = []


    for event_id in event_dict:
        hits = event_dict[event_id]
        no = int(event_id) + 1
        if no in event_dict:
            next_hit = event_dict[no]
            next_exists = True
        else:
            next_exists = False
        
        if len(hits) == 1:
            if next_exists and find_true((hits[0], next_hit[0])):
                group.append(((hits[0], next_hit[0]), "true"))
                #LOR_list.append(find_LOR((hits[0], next_hit[0])))
            else:
                group.append((hits[0], "random"))
        elif len(hits) == 2:
            e1,e2 = hits[0], hits[1]
            
            if find_true((e1,e2)):
                group.append(((e1,e2), "true"))
                LOR_list.append(find_LOR((e1,e2)))
            else:
                group.append(((e1,e2), "scatter"))
        else:
                
                complex_no += len(hits)
                group.append((hits, "complex"))

        
    return group, complex_no, LOR_list

def get_data(file):
    """
    finds necr data (and true,scatter,random, complex counts). estimates decay site
    """

    grouped_results, complex_counts, LORs = group_counts(file)
    detected_counts =len(file)

    true_counts = sum(2 for result in grouped_results if result[1] == "true") #each event counts for 2
    scatter_counts = sum(2 for result in grouped_results if result[1] == "scatter") # same here
    random_counts = sum(1 for result in grouped_results if result[1] == "random") #single randoms

    midpoints = [find_midpoint(lor) for lor in LORs]
    event_coords = find_decay_site(midpoints)

    necr = find_necr(true_counts, scatter_counts, random_counts, 1)

    return detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords

def file_reader(file_path):
    """
    makes csv file into numpy array for analysis
    """
    df = pd.read_csv(file_path).to_numpy()
    return df

def analyse_data():
    """
    plots/prints all data analysis results
    """

    necr_list = []
    trues_list = []
    scatters_list = []
    randoms_list = []
    detected_list = []
    complex_list = []
    event_coords_list = []

    necr_rates = []
    dose = 400 #MBq

    #loop over all files and collect data
    for i in range(27):
        file_path = "data/" + files[i] + ".csv"
        data = file_reader(file_path)
        detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
        complex_list.append(complex_counts)
        necr_list.append(necr)
        trues_list.append(true_counts)
        scatters_list.append(scatter_counts)
        randoms_list.append(random_counts)
        detected_list.append(detected_counts)
        event_coords_list.append(event_coords)

        necr_rate = necr_as_rate(necr, dose, detected_counts)
        necr_rates.append(necr_rate)
    
    
    #for estimated decay site analysis:
    xi_list, yi_list, zi_list = zip(*event_coords_list)

    dxi_list = [abs(xi - actual_site[0]) for xi in xi_list]
    dyi_list = [abs(yi - actual_site[1]) for yi in yi_list]
    dzi_list = [abs(zi - actual_site[2]) for zi in zi_list]


    #print site with least difference from actual x coord... can be repeated for y and z if needed
    #min_diff_index = dxi_list.index(min(dxi_list))
    #print(f"Detector length with least x difference from actual site: {detector_lengths_mm[min_diff_index]} mm, Difference: {min(dxi_list)} mm")
    
    #average coordinates across all lengths... as they are all kind of the same anyways
    avg_xi = sum(xi_list) / len(xi_list)
    avg_yi = sum(yi_list) / len(yi_list)
    avg_zi = sum(zi_list) / len(zi_list)

    #uncertainties on averages.. ish
    uncert_x = max(xi_list) - min(xi_list)
    uncert_y = max(yi_list) - min(yi_list)
    uncert_z = max(zi_list) - min(zi_list)

    print(f"Average estimated decay site: (x={avg_xi:.2f} ± {uncert_x/2:.2f}, y={avg_yi:.2f} ± {uncert_y/2:.2f}, z={avg_zi:.2f} ± {uncert_z/2:.2f})")          
    
    #finds difference from actual site for each estimate... hashed out cause it clutters output
    #for i in range(len(detector_lengths_mm)):
        #xi, yi, zi = event_coords_list[i]

        #print(f"Length {detector_lengths_mm[i]} mm: (x={xi:.2f}, y={yi:.2f}, z={zi:.2f})")
        #dxi, dyi, dzi = compute_decay_site_difference(event_coords_list[i], actual_site)
        #print(f"Difference from actual site: (dx={dxi:.2f}, dy={dyi:.2f}, dz={dzi:.2f})")
        #print(f"Difference from actual site: {compute_decay_site_difference(event_coords_list[i], actual_site)}")

    plot_necr(necr_list, detector_lengths_mm)
    plot_necr_rates(necr_rates, detector_lengths_mm)
    plot_trues(trues_list, detector_lengths_mm)
    plot_all_counts(trues_list, scatters_list, randoms_list, complex_list, detector_lengths_mm)
    plot_total_counts(detected_list, detector_lengths_mm)
    plot_complex(complex_list, detector_lengths_mm)
    plot_site_coords(xi_list, yi_list, zi_list, detector_lengths_mm, actual_site)


def full_body_decays():
    filepath = "data/wholebodyhits.csv"
    data = file_reader(filepath)
    detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
    total_counts = detected_counts + true_counts + scatter_counts + random_counts + complex_counts
    
    print("Total detected counts from full body source data:", total_counts)
    print("Estimated decay site from full body source data:", event_coords)

    #print(f"Detected counts from full body source {detected_counts}")
    #print(f"True counts from full body source {true_counts}")
    #print(f"Scatter counts from full body source {scatter_counts}")
    #print(f"Random counts from full body source {random_counts}")
    #print(f"Complex counts from full body source {complex_counts}")
    #print(f"NECR from full body source {necr}")

def z_axis_decays():
    """
    note: this data is from the cylindrical phantom so NECR values will not be comparable to other data, 
    only rly use this for decay site comparison
    """
    filepath = "data/hits_linear_source.csv"
    data = file_reader(filepath)
    detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
    total_counts = detected_counts + true_counts + scatter_counts + random_counts + complex_counts

    print("Total detected counts from z-axis source data:", total_counts)
    print("Estimated decay site from z-axis source data:", event_coords)

def both_breast_decays():
    filepath = "data/bothbreasthits.csv"
    data = file_reader(filepath)
    detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
    total_counts = detected_counts + true_counts + scatter_counts + random_counts + complex_counts

    print("Total detected counts from both breast source data:", total_counts)
    print("Estimated decay site from both breast source data:", event_coords)

    #print(f"Detected counts from both breast source {detected_counts}")
    #print(f"True counts from both breast source {true_counts}")
    #print(f"Scatter counts from both breast source {scatter_counts}")
    #print(f"Random counts from both breast source {random_counts}")
    #print(f"Complex counts from both breast source {complex_counts}")
    #print(f"NECR from both breast source {necr}")

def control_data():
    """
    run control data analyses
    """
    full_body_decays()
    z_axis_decays()
    both_breast_decays()

# run all our analyses
analyse_data()
control_data()