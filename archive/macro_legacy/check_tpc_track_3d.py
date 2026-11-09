#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Macro to check TPC tracking results in 3D: interactive 3D visualization using matplotlib
Usage:
  python check_tpc_track_3d.py <rootfile>
  or
  python
  >>> from check_tpc_track_3d import *
  >>> set_path_3d("path/to/rootfile.root")
  >>> event_3d()        # random event
  >>> event_3d(evnum)   # specific event

Controls:
  - Left mouse button + drag: Rotate view
  - Right mouse button + drag: Zoom
  - Middle mouse button + drag: Pan
  - Mouse wheel: Zoom in/out
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import ROOT
import random

# Global variables
g_file_3d = None
g_tree_3d = None
g_entries_3d = 0
g_fig_3d = None
g_ax_3d = None

# Event data structure
class EventData3D:
    def __init__(self):
        self.runnum = 0
        self.evnum = 0
        self.nhTpc = 0
        self.nclTpc = 0
        self.ntTpc = 0
        self.raw_hitpos_x = []
        self.raw_hitpos_y = []
        self.raw_hitpos_z = []
        self.raw_de = []
        self.raw_padid = []
        self.raw_layer = []
        self.raw_row = []
        self.cluster_x = []
        self.cluster_y = []
        self.cluster_z = []
        self.cluster_de = []
        self.cluster_layer = []
        self.cluster_row_center = []
        self.cluster_houghflag = []
        self.x0Tpc = []
        self.y0Tpc = []
        self.u0Tpc = []
        self.v0Tpc = []
        self.theta = []
        self.nhtrack = []
        self.hitpos_x = []
        self.hitpos_y = []
        self.hitpos_z = []

g_event_3d = EventData3D()

def set_path_3d(path):
    """Open ROOT file and prepare tree for reading"""
    global g_file_3d, g_tree_3d, g_entries_3d
    
    if g_file_3d:
        g_file_3d.Close()
    
    g_file_3d = ROOT.TFile.Open(path)
    if not g_file_3d or g_file_3d.IsZombie():
        print(f"Error: Cannot open file {path}")
        return
    
    g_tree_3d = g_file_3d.Get("tpc")
    if not g_tree_3d:
        print("Error: Cannot find tree 'tpc' in file")
        return
    
    g_entries_3d = g_tree_3d.GetEntries()
    print(f"File opened: {path}")
    print(f"Total entries: {g_entries_3d}")

def load_event_3d(entry=-1):
    """Load event data from ROOT tree"""
    global g_tree_3d, g_entries_3d, g_event_3d
    
    if not g_tree_3d:
        print("Error: No file opened. Use set_path_3d() first.")
        return
    
    if g_entries_3d == 0:
        print("Error: No entries in tree")
        return
    
    if entry < 0:
        entry = random.randint(0, g_entries_3d - 1)
    
    if entry >= g_entries_3d:
        print(f"Error: Entry {entry} is out of range [0, {g_entries_3d})")
        return
    
    g_tree_3d.GetEntry(entry)
    
    # Get branch values
    g_event_3d.runnum = g_tree_3d.run_number
    g_event_3d.evnum = g_tree_3d.event_number
    g_event_3d.nhTpc = g_tree_3d.nhTpc
    g_event_3d.nclTpc = g_tree_3d.nclTpc
    g_event_3d.ntTpc = g_tree_3d.ntTpc
    
    # Convert ROOT vectors to Python lists
    g_event_3d.raw_hitpos_x = list(g_tree_3d.raw_hitpos_x)
    g_event_3d.raw_hitpos_y = list(g_tree_3d.raw_hitpos_y)
    g_event_3d.raw_hitpos_z = list(g_tree_3d.raw_hitpos_z)
    g_event_3d.raw_de = list(g_tree_3d.raw_de)
    g_event_3d.raw_padid = list(g_tree_3d.raw_padid)
    g_event_3d.raw_layer = list(g_tree_3d.raw_layer)
    g_event_3d.raw_row = list(g_tree_3d.raw_row)
    
    g_event_3d.cluster_x = list(g_tree_3d.cluster_x)
    g_event_3d.cluster_y = list(g_tree_3d.cluster_y)
    g_event_3d.cluster_z = list(g_tree_3d.cluster_z)
    g_event_3d.cluster_de = list(g_tree_3d.cluster_de)
    g_event_3d.cluster_layer = list(g_tree_3d.cluster_layer)
    g_event_3d.cluster_row_center = list(g_tree_3d.cluster_row_center)
    g_event_3d.cluster_houghflag = list(g_tree_3d.cluster_houghflag)
    
    g_event_3d.x0Tpc = list(g_tree_3d.x0Tpc)
    g_event_3d.y0Tpc = list(g_tree_3d.y0Tpc)
    g_event_3d.u0Tpc = list(g_tree_3d.u0Tpc)
    g_event_3d.v0Tpc = list(g_tree_3d.v0Tpc)
    g_event_3d.theta = list(g_tree_3d.theta)
    g_event_3d.nhtrack = list(g_tree_3d.nhtrack)
    
    # Handle nested vectors
    g_event_3d.hitpos_x = []
    g_event_3d.hitpos_y = []
    g_event_3d.hitpos_z = []
    for i in range(g_event_3d.ntTpc):
        if i < len(g_tree_3d.hitpos_x):
            g_event_3d.hitpos_x.append(list(g_tree_3d.hitpos_x[i]))
            g_event_3d.hitpos_y.append(list(g_tree_3d.hitpos_y[i]))
            g_event_3d.hitpos_z.append(list(g_tree_3d.hitpos_z[i]))
        else:
            g_event_3d.hitpos_x.append([])
            g_event_3d.hitpos_y.append([])
            g_event_3d.hitpos_z.append([])
    
    print(f"Event loaded: Run={g_event_3d.runnum}, Event={g_event_3d.evnum}, Entry={entry}")
    print(f"  Raw hits: {g_event_3d.nhTpc}")
    print(f"  Clusters: {g_event_3d.nclTpc}")
    print(f"  Tracks: {g_event_3d.ntTpc}")

def event_3d(evnum=-1):
    """Display 3D visualization of TPC tracks"""
    global g_fig_3d, g_ax_3d, g_event_3d
    
    load_event_3d(evnum)
    
    # Create or clear figure
    if g_fig_3d is None:
        g_fig_3d = plt.figure(figsize=(12, 8))
        g_ax_3d = g_fig_3d.add_subplot(111, projection='3d')
    else:
        g_ax_3d.clear()
    
    # Calculate range
    xmin, xmax = -300.0, 300.0
    ymin, ymax = -200.0, 200.0
    zmin, zmax = -300.0, 300.0
    
    # Find actual range from data
    if g_event_3d.nhTpc > 0:
        if len(g_event_3d.raw_hitpos_x) > 0:
            xmin = min(g_event_3d.raw_hitpos_x)
            xmax = max(g_event_3d.raw_hitpos_x)
        if len(g_event_3d.raw_hitpos_y) > 0:
            ymin = min(g_event_3d.raw_hitpos_y)
            ymax = max(g_event_3d.raw_hitpos_y)
        if len(g_event_3d.raw_hitpos_z) > 0:
            zmin = min(g_event_3d.raw_hitpos_z)
            zmax = max(g_event_3d.raw_hitpos_z)
    
    # Add margin
    xmargin = (xmax - xmin) * 0.1
    ymargin = (ymax - ymin) * 0.1
    zmargin = (zmax - zmin) * 0.1
    xmin -= xmargin; xmax += xmargin
    ymin -= ymargin; ymax += ymargin
    zmin -= zmargin; zmax += zmargin
    
    # Draw raw hits
    if g_event_3d.nhTpc > 0 and len(g_event_3d.raw_hitpos_x) > 0:
        g_ax_3d.scatter(g_event_3d.raw_hitpos_x,
                       g_event_3d.raw_hitpos_y,
                       g_event_3d.raw_hitpos_z,
                       c='gray', marker='o', s=10, alpha=0.5, label='Raw hits')
    
    # Draw clusters
    if g_event_3d.nclTpc > 0 and len(g_event_3d.cluster_x) > 0:
        g_ax_3d.scatter(g_event_3d.cluster_x,
                       g_event_3d.cluster_y,
                       g_event_3d.cluster_z,
                       c='blue', marker='^', s=50, alpha=0.7, label='Clusters')
    
    # Draw Hough-selected clusters
    if g_event_3d.nclTpc > 0 and len(g_event_3d.cluster_houghflag) > 0:
        hough_x = []
        hough_y = []
        hough_z = []
        for i in range(min(g_event_3d.nclTpc, len(g_event_3d.cluster_houghflag))):
            if g_event_3d.cluster_houghflag[i] > 0:
                if i < len(g_event_3d.cluster_x):
                    hough_x.append(g_event_3d.cluster_x[i])
                    hough_y.append(g_event_3d.cluster_y[i])
                    hough_z.append(g_event_3d.cluster_z[i])
        if len(hough_x) > 0:
            g_ax_3d.scatter(hough_x, hough_y, hough_z,
                           c='red', marker='s', s=80, alpha=0.8, label='Hough selected')
    
    # Draw tracks
    colors = plt.cm.tab10(np.linspace(0, 1, max(g_event_3d.ntTpc, 1)))
    for itrack in range(g_event_3d.ntTpc):
        if (itrack < len(g_event_3d.x0Tpc) and
            itrack < len(g_event_3d.y0Tpc) and
            itrack < len(g_event_3d.u0Tpc) and
            itrack < len(g_event_3d.v0Tpc)):
            x0 = g_event_3d.x0Tpc[itrack]
            y0 = g_event_3d.y0Tpc[itrack]
            u0 = g_event_3d.u0Tpc[itrack]
            v0 = g_event_3d.v0Tpc[itrack]
            
            # Calculate track points at z boundaries
            z1 = zmin
            z2 = zmax
            x1 = x0 + u0 * z1
            y1 = y0 + v0 * z1
            x2 = x0 + u0 * z2
            y2 = y0 + v0 * z2
            
            label = f'Track {itrack}' if itrack == 0 else None
            g_ax_3d.plot([x1, x2], [y1, y2], [z1, z2],
                        color=colors[itrack % len(colors)], linewidth=2, label=label)
    
    # Draw track hits if available
    for itrack in range(g_event_3d.ntTpc):
        if (itrack < len(g_event_3d.nhtrack) and
            g_event_3d.nhtrack[itrack] > 0 and
            itrack < len(g_event_3d.hitpos_x)):
            nhits = g_event_3d.nhtrack[itrack]
            if (nhits > 0 and
                nhits <= len(g_event_3d.hitpos_x[itrack])):
                g_ax_3d.scatter(g_event_3d.hitpos_x[itrack][:nhits],
                               g_event_3d.hitpos_y[itrack][:nhits],
                               g_event_3d.hitpos_z[itrack][:nhits],
                               c=[colors[itrack % len(colors)]], marker='*', s=100, alpha=0.9)
    
    # Set labels and limits
    g_ax_3d.set_xlabel('X [mm]')
    g_ax_3d.set_ylabel('Y [mm]')
    g_ax_3d.set_zlabel('Z [mm]')
    g_ax_3d.set_xlim(xmin, xmax)
    g_ax_3d.set_ylim(ymin, ymax)
    g_ax_3d.set_zlim(zmin, zmax)
    g_ax_3d.set_title(f'TPC Track 3D View (Run {g_event_3d.runnum}, Event {g_event_3d.evnum})')
    g_ax_3d.legend()
    
    plt.tight_layout()
    plt.show()
    
    # Print summary
    print("\n=== Event Summary (3D) ===")
    print(f"Raw hits: {g_event_3d.nhTpc}")
    print(f"Clusters: {g_event_3d.nclTpc}")
    nHoughClusters = sum(1 for f in g_event_3d.cluster_houghflag if f > 0)
    print(f"Hough-selected clusters: {nHoughClusters}")
    print(f"Tracks: {g_event_3d.ntTpc}")
    for itrack in range(g_event_3d.ntTpc):
        if itrack < len(g_event_3d.x0Tpc):
            nhits = g_event_3d.nhtrack[itrack] if itrack < len(g_event_3d.nhtrack) else 0
            theta_val = g_event_3d.theta[itrack] if itrack < len(g_event_3d.theta) else 0.0
            print(f"  Track {itrack}: x0={g_event_3d.x0Tpc[itrack]:.2f}, "
                  f"y0={g_event_3d.y0Tpc[itrack]:.2f}, "
                  f"u0={g_event_3d.u0Tpc[itrack]:.4f}, "
                  f"v0={g_event_3d.v0Tpc[itrack]:.4f}, "
                  f"theta={theta_val:.2f}, nhits={nhits}")
    print("=========================")
    print("\n3D View Controls:")
    print("  - Left mouse + drag: Rotate")
    print("  - Right mouse + drag: Zoom")
    print("  - Middle mouse + drag: Pan")
    print("  - Mouse wheel: Zoom in/out")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_tpc_track_3d.py <rootfile> [entry]")
        sys.exit(1)
    
    set_path_3d(sys.argv[1])
    entry = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    event_3d(entry)
