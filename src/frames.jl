function align_frame_with_vector(target_primary, target_secondary, axis_primary, axis_secondary)
    target_primary = normalize(target_primary)
    target_secondary = normalize(target_secondary)
    axis1 = normalize(cross(axis_primary, target_primary))
    angle1 = acos(dot(axis_primary, target_primary))
    q1 = Quaternion(axis1, angle1)

    secondary_perpendicular = normalize(target_secondary - dot(target_secondary, target_primary) * target_primary)
    secondary_intermediate = rotvec(axis_secondary, q1)
    axis2 = normalize(cross(secondary_intermediate, secondary_perpendicular))
    angle2 = acos(dot(secondary_intermediate, secondary_perpendicular))
    q2 = Quaternion(axis2, angle2)
    
    q = q2 * q1
end 