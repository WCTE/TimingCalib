# -*- coding: utf-8 -*-
"""
Device: A base class for parts that make up a water Cherenkov detector system

@author: karlen
"""
import numpy as np
from scipy import stats
from scipy.spatial.transform import Rotation as R

class Device:
    """A base class for parts that make up a detector system
    The following are accessible:
        - property dictionaries
            .prop_design: intended properties (read-only)
            .prop_true: actual properties (read-only)
            .prop_est: current property estimates
            .prop_est_sig: standard deviation of property estimators
            .prop_prior: prior property estimates (TBD)
            .prop_prior_sig: standard deviation of priors (TBD)
        - placement dictionaries:
            .place_design: intended placements (read-only)
            .place_true: actual placement (read-only)
            .place_est: current placement estimates
            .place_est_sig: standard deviation of placement estimators
    Misc:
        - "device_type" points to a device class (e.g. PMT or LED)
        - "kind" defines the particular version of the device_type
        - "type" is a reserved word in python - avoid using this as a variable name
        - when building a device of a particular kind, the device class
        properties "prop_mean, prop_scale, prop_var" define the actual properties
        of the device. These class properties are dictionaries keyed by "kind"
        - "prop_mean" is the mean value
        - "prop_scale" is the scale of the variation (sigma for normal,
        half-width for uniform). If 0, there is no variation
        - "prop_var" specifies variation type if not Normal. "norm" or "uniform"
        - "rot_axes" is the combination of up to 3 Euler rotations, e.g. 'XZX'
        look at scipy.spatial.rotation documentation for conventions
        - "rot_angles" is the corresponding rotation angles in radians
        - "rot_angles_sig" is the corresponding standard deviations of the angles

    """

    def __init__(self):
        """Constructor"""

        self._prop_design = {}
        self._prop_true = {}
        self.prop_est = {}
        self.prop_est_sig = {}

        self._place_design = {}
        self._place_true = {}
        self.place_est = {}
        self.place_est_sig = {}

        self.kind = None
        self.controller = None

    def set_properties(self, device_type, kind):
        """Set the instance properties for the device object"""
        if kind in device_type.prop_mean:
            self._prop_design = device_type.prop_mean[kind]
            truth = {}
            for key in device_type.prop_mean[kind]:
                scale = (device_type.prop_scale[kind])[key]
                mean = (device_type.prop_mean[kind])[key]
                if scale > 0.:
                    var_type = (device_type.prop_var[kind]).get(key, 'norm')
                    if var_type == 'uniform':
                        val = stats.uniform.rvs(loc=mean-scale, scale=2.*scale)
                    else:
                        val = stats.norm.rvs(loc=mean, scale=scale)
                    truth[key] = val
                else:
                    truth[key] = mean
            self._prop_true = truth

    def set_placement(self, place_design, place_true):
        """Set the device placement information"""
        self._place_design = place_design.copy()
        self._place_true = place_true.copy()

    def place_devices(self, device_type, devices_design, kind):
        """Places the sub_devices that make up this device"""
        devices = []
        if kind in devices_design:
            for device in devices_design[kind]:
                device_kind = device['kind']
                place_design = {'loc':device['loc'],
                                'rot_axes':device['rot_axes'],
                                'rot_angles':device['rot_angles']}
                loc_true = []
                for i in range(3):
                    loc_true.append(stats.norm.rvs(device['loc'][i], device['loc_sig'][i]))

                rot_angles_true = []
                for i in range(len(device['rot_angles'])):
                    rot_angles_true.append(stats.norm.rvs(device['rot_angles'][i],
                                                          device['rot_angles_sig'][i]))
                place_true = {'loc':loc_true, 'rot_axes':device['rot_axes'],
                              'rot_angles':rot_angles_true}
                new_device = device_type(self, device_kind, place_design, place_true)
                devices.append(new_device)
        return devices

    def get_z_orientation(self, place_info):
        """Return the location af a device and direction of its z axis in WCD coordinate system"""
        # place_info is a string: either 'true', 'design', or 'est'
        # head is the location of the device, tail is the end of a vector of
        # length dist_head aligned with the device's z-axis
        device_place = getattr(self, 'place_'+place_info, None)
        head = device_place['loc']
        dist_head = 100.

        rot1 = R.from_euler(device_place['rot_axes'], device_place['rot_angles'])
        z_axis = rot1.apply([0., 0., 1.])

        if self.__class__.__name__ == 'MPMT':
            location = head
            direction = z_axis
        else:
            # apply translation and rotation of mPMT
            tail = np.add(head, dist_head*z_axis)
            mpmt_place = getattr(self.controller, 'place_'+place_info, None)
            rot_head = head
            rot_tail = tail
            if 'rot_axes' in mpmt_place:
                rot2 = R.from_euler(mpmt_place['rot_axes'],
                                    mpmt_place['rot_angles'])
                rot_head = rot2.apply(head)
                rot_tail = rot2.apply(tail)
            new_direction = np.subtract(rot_tail, rot_head)
            norm = np.linalg.norm(new_direction)

            location = np.add(rot_head, mpmt_place.get('loc',[0.,0.,0.]))
            direction = np.divide(new_direction, norm)

        return location, direction

    def get_orientation(self, place_info):
        """Return the location af device and directions of its x and z axes in WCD coordinate system"""
        # place_info is a string: either 'true', 'design', or 'est'
        # head is the location of the device, tail_x,_z are the ends of vectors of
        # length dist_head aligned with the device's x and z-axes
        device_place = getattr(self, 'place_'+place_info, None)
        head = device_place['loc']
        dist_head = 100.

        rot1 = R.from_euler(device_place['rot_axes'], device_place['rot_angles'])
        x_axis = rot1.apply([1., 0., 0.])
        z_axis = rot1.apply([0., 0., 1.])

        if self.__class__.__name__ == 'MPMT':
            location = head
            direction_x = x_axis
            direction_z = z_axis
        else:
            # apply translation and rotation of mPMT
            tail_x = np.add(head, dist_head * x_axis)
            tail_z = np.add(head, dist_head * z_axis)
            mpmt_place = getattr(self.controller, 'place_'+place_info, None)
            rot_head = head
            rot_tail_x = tail_x
            rot_tail_z = tail_z
            if 'rot_axes' in mpmt_place:
                rot2 = R.from_euler(mpmt_place['rot_axes'], mpmt_place['rot_angles'])
                rot_head = rot2.apply(head)
                rot_tail_x = rot2.apply(tail_x)
                rot_tail_z = rot2.apply(tail_z)
            new_direction_x = np.subtract(rot_tail_x, rot_head)
            new_direction_z = np.subtract(rot_tail_z, rot_head)
            norm_x = np.linalg.norm(new_direction_x)
            norm_z = np.linalg.norm(new_direction_z)

            location = np.add(rot_head, mpmt_place.get('loc',[0.,0.,0.]))
            direction_x = np.divide(new_direction_x, norm_x)
            direction_z = np.divide(new_direction_z, norm_z)

        return location, direction_x, direction_z

    def get_circle_points(self, n_point, place_info):
        """Return a list of length n_point space points of the circle in the x-y plane"""
        device_prop = getattr(self, 'prop_' + place_info, None)
        radius = device_prop['size']/2.
        location, direction_x, direction_z = self.get_orientation(place_info)

        # change magnitude of direction_x to device radius
        perp = np.array(direction_x)*radius

        # small rotation about the device z_axis
        rot = R.from_rotvec(2. * np.pi / n_point * np.array(direction_z))
        points = []
        for i in range(n_point):
            points.append(list(np.add(location, perp)))
            perp = rot.apply(perp)

        return points


    @property
    def prop_true(self):
        """access method for read-only data"""
        return self._prop_true

    @property
    def prop_design(self):
        """access method for read-only data"""
        return self._prop_design

    @property
    def place_true(self):
        """access method for read-only data"""
        return self._place_true

    @property
    def place_design(self):
        """access method for read-only data"""
        return self._place_design
