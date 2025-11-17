import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


EVENT_ID, GEOMETRY_ID, ENERGY, TIME, R, PHI, Z = 0, 1, 2, 3, 4, 5, 6 #units: keV, ns, mm, rad
detector_lengths_mm = [100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050, 1100,1150,1200,1250,1300,1350,1400]
files = ["hits100mm", "hits150mm", "hits200mm", "hits250mm", "hits300mm", \
"hits350mm", "hits400mm", "hits450mm", "hits500mm", "hits550mm", "hits600mm", "hits650mm", "hits700mm", "hits750mm", \
"hits800mm", "hits850mm", "hits900mm", "hits950mm", "hits1000mm", "hits1050mm", "hits1100mm", "hits1150mm", "hits1200mm", \
"hits1250mm", "hits1300mm", "hits1350mm", "hits1400mm"]

mm100, mm150, mm200, mm250, mm300, mm350, mm400, mm450, mm500, mm550, mm600, mm650, mm700, \
mm750, mm800, mm850, mm900, mm950, mm1000, mm1050 = "hits100mm", "hits150mm", "hits200mm", "hits250mm", "hits300mm", \
"hits350mm", "hits400mm", "hits450mm", "hits500mm", "hits550mm", "hits600mm", "hits650mm", "hits700mm", "hits750mm", \
"hits800mm", "hits850mm", "hits900mm", "hits950mm", "hits1000mm", "hits1050mm"

simulated_counts = 100000  #number of simulated decays per run

def find_true(event):
    """Determine if an event is a true coincidence.
    Takes a touple of two detections (event) and checks if both energies are within 50kev of 511kev.
    Checks whether photons are unscattered i.e. angle = 180 degrees (within 5 degrees).
    """
    energy1 = float(event[0][ENERGY])
    energy2 = float(event[1][ENERGY])
    angle = find_angle(event)
    
    
    is_true = (abs(energy1 - 511) < 70) and (abs(energy2 - 511) < 70) #and (abs(abs(angle) - 180) < 50)
    return is_true

def find_scatter(e1,e2):
    
    energy1 = float(e1[ENERGY])
    energy2 = float(e2[ENERGY])
    angle = find_angle((e1,e2))
    
    
    
    return is_scatter

def find_angle(event):
    """Calculate the angle between two hits."""
    phi0 = float(event[0][PHI])
    phi1 = float(event[1][PHI])
    r0 = float(event[0][R])
    r1 = float(event[1][R])
    #print("type of phi0 = ", type(phi0))

    x1 = r0*np.cos(phi0)
    y1 = r0*np.sin(phi0)
    
    x2 = r1*np.cos(phi1)
    y2 = r1*np.sin(phi1)
    
    delta_x = x2 - x1
    delta_y = y2 - y1
    
    angle = np.arctan2(delta_y, delta_x) * (180 / np.pi)  # Convert to degrees
    #print("angle = ", angle)
    return angle


def determine_four_complex(events):
    "Return number of true coincidences found in all combinations of two hits in events"
    true_count = 0
    for i in range(len(events)):
        for j in range(i+1, len(events)):
            if find_true((events[i], events[j])):
                true_count += 1
    return true_count
def determine_three_complex(events):
    "Return True if any combination of two hits in events is a true coincidence. False otherwise"
    for i in range(len(events)):
        for j in range(i+1, len(events)):
            if find_true((events[i], events[j])):
                return True
    return False


def group_counts(df):
    """Group counts by event ID, and determine whether each event is true, scatter, or random."""
    #find all events with same event id
    #if 1 - random 
    #if 2 - check if true or scatter
    #if >2 - random for now.. 

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
        #print("event_id = ", event_id)
        hits = event_dict[event_id]
        no = int(event_id) + 1
        if no in event_dict:
            next_hit = event_dict[no]
            next_exists = True
        else:
            next_exists = False
        
        #print("number of hits = ", len(hits))
        if len(hits) == 1:
            if next_exists and find_true((hits[0], next_hit[0])):
                group.append(((hits[0], next_hit[0]), "true"))
                #LOR_list.append(find_LOR((hits[0], next_hit[0])))
            else:
                group.append((hits[0], "random"))
        elif len(hits) == 2:
            e1,e2 = hits[0], hits[1]
            #print(f" hits1= ", e1)
            #print(f" hits2= ", e2)
            
            if find_true((e1,e2)):
                group.append(((e1,e2), "true"))
                LOR_list.append(find_LOR((e1,e2)))
            else:
                group.append(((e1,e2), "scatter"))
        else:
                
                complex_no += len(hits)
                group.append((hits, "scatter"))
                
                    
                   # print("Event ID with more than 3 hits detected:", event_id)
                """
            for i in range(len(hits)):
                for j in range(i+1, len(hits)):
                    if find_true((hits[i], hits[j])):
                        group.append(((hits[i], hits[j]), "true"))

                    else:
                        group.append(((hits[i], hits[j]), "scatter"))
                        #need to work this one out properly
            #group.append((hits, "random"))
            #check if 2 of detected counts are true/scatter
            """

        
    return group, complex_no, LOR_list        
        
def find_necr(t,s,r,k):
    """Calculate NECR given true (t), scatter (s), random (r) counts, and scale factor k."""
    return (t*t) / (t + s + k*r)

def get_data(file):
    
    grouped_results, complex_counts, LORs = group_counts(file)
    detected_counts =len(file)
    true_counts = sum(2 for result in grouped_results if result[1] == "true") #each event counts for 2.. or 1 cause both events added?
    scatter_counts = sum(2 for result in grouped_results if result[1] == "scatter") # same here
    random_counts = sum(1 for result in grouped_results if result[1] == "random") #single randoms, definitely one
    midpoints = [find_midpoint(lor) for lor in LORs]
    event_coords = find_decay_site(midpoints)
    #complex_counts = sum(3 for result in grouped_results if result[1] == "complexcase")
    print(f"Sum of grouped results: {true_counts + scatter_counts + random_counts + complex_counts}")
    print(f"Complex counts (more than 2 hits per event): {complex_counts}")
    necr = find_necr(true_counts, scatter_counts, random_counts, 1)
    return detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords

def file_reader(file_path):
    """Read data from a file and return as a list of lines."""
    df = pd.read_csv(file_path).to_numpy()
    return df

def compare_data():
    necr_list = []
    trues_list = []
    scatters_list = []
    randoms_list = []
    detected_list = []
    complex_list = []
    event_coords_list = []

    for i in range(27):
        file_path = "/Users/hedda/Desktop/Y4/SH Project/CODE/csv files/" + files[i] + ".csv"
        data = file_reader(file_path)
        detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
        complex_list.append(complex_counts)
        necr_list.append(necr)
        trues_list.append(true_counts)
        scatters_list.append(scatter_counts)
        randoms_list.append(random_counts)
        detected_list.append(detected_counts)
        event_coords_list.append(event_coords)
    
    print("Event coordinates estimates for each detector length:")
    for i in range(len(detector_lengths_mm)):
        xi, yi, zi = event_coords_list[i]

        print(f"Length {detector_lengths_mm[i]} mm: (x={xi:.2f}, y={yi:.2f}, z={zi:.2f})")
        dxi, dyi, dzi = compute_decay_site_difference(event_coords_list[i], actual_site)
        print(f"Difference from actual site: (dx={dxi:.2f}, dy={dyi:.2f}, dz={dzi:.2f})")
        #print(f"Difference from actual site: {compute_decay_site_difference(event_coords_list[i], actual_site)}")

    plot_necr(necr_list, detector_lengths_mm)
    plot_trues(trues_list, detector_lengths_mm)
    plot_all_counts(trues_list, scatters_list, randoms_list, complex_list, detector_lengths_mm)
    plot_total_counts(detected_list, detector_lengths_mm)


def plot_necr(ns, lengths):
    plt.plot(lengths, ns, marker='o')
    plt.title('NECR vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('NECR')
    
    plt.show()

def plot_trues(ts, lengths):
    plt.plot(lengths, ts, marker='o')
    plt.title('True Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('True Counts')
    
    plt.show()

def plot_all_counts(ts,ss,rs,cs,lengths):
    plt.plot(lengths, ts, marker='o', label='True Counts')
    plt.plot(lengths, ss, marker='o', label='Scatter Counts')
    plt.plot(lengths, rs, marker='o', label='Single Counts')
    plt.plot(lengths, cs, marker='o', label='Complex (>2 hits) Counts')
    plt.title('Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Counts')
    plt.legend()
    
    plt.show()
    
def plot_total_counts(totals,lengths):
    #plot straight line of simulated counts for reference
    plt.axhline(y=simulated_counts, color='r', linestyle='--', label='Simulated Decays')
    plt.plot(lengths, totals, marker='o', label='Total Detected Counts')
    plt.title('Total Detected Counts vs Detector Length')
    plt.xlabel('Detector Length (mm)')
    plt.ylabel('Total Detected Counts')
    plt.legend()
    
    plt.show()

def find_LOR(event):
    "calculate line of response from two hits"
    phi0 = float(event[0][PHI])
    phi1 = float(event[1][PHI])
    r0 = float(event[0][R])
    r1 = float(event[1][R])
    #print("type of phi0 = ", type(phi0))

    x1 = r0*np.cos(phi0)
    y1 = r0*np.sin(phi0)
    
    x2 = r1*np.cos(phi1)
    y2 = r1*np.sin(phi1)

    z1 = float(event[0][Z])
    z2 = float(event[1][Z])

    return (x1, y1, z1), (x2, y2, z2)

def find_midpoint(line):
    "find midpoint of a line defined by two endpoints"
    (x1, y1, z1), (x2, y2, z2) = line
    x_mid = (x1 + x2) / 2
    y_mid = (y1 + y2) / 2
    z_mid = (z1 + z2) / 2
    return (x_mid, y_mid, z_mid)

def find_decay_site(midpoints):
    n = len(midpoints)
    x_sum = 0
    y_sum = 0
    z_sum = 0
    for midpoint in midpoints:
        x_sum += midpoint[0]
        y_sum += midpoint[1]
        z_sum += midpoint[2]
    return (x_sum / n, y_sum / n, z_sum / n)

actual_site = (200.3/4, 100.0/2, 0.0)

def compute_decay_site_difference(estimated, actual):
    "compute difference between estimated and actual decay site"
    dx = abs(estimated[0] - actual[0])
    dy = abs(estimated[1] - actual[1])
    dz = abs(estimated[2] - actual[2])
    return (dx, dy, dz)   

def linear_source_data():
    filepath = "/Users/hedda/Desktop/Y4/SH Project/CODE/csv files/hits_linear_source.csv"
    data = file_reader(filepath)
    detected_counts, true_counts, scatter_counts, random_counts, complex_counts, necr, event_coords = get_data(data)
    
    
    
    print("Estimated decay site from linear source data:", event_coords)

    print(f"Detected counts from linear source, cylindrical phantom: {detected_counts}")
    print(f"True counts from linear source, cylindrical phantom: {true_counts}")
    print(f"Scatter counts from linear source, cylindrical phantom: {scatter_counts}")
    print(f"Random counts from linear source, cylindrical phantom: {random_counts}")
    print(f"Complex counts from linear source, cylindrical phantom: {complex_counts}")
    print(f"NECR from linear source, cylindrical phantom: {necr}")

compare_data()
linear_source_data()


test = "1935 903 510.999 1.73967 419.356 -1.46742 222.855"
test2 = "1935 957 510.999 1.47194 413.229 1.08676 279.187"


t = find_true((test.split(), test2.split()))
print("is true: ", t)