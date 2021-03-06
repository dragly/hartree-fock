#version 330 core
uniform sampler3D myTexture3D;
uniform vec4 ve_eyePosition; // last column of the inverse modelViewMatrix
uniform float multiplier;
uniform float quality;
uniform float contrast;
uniform float mixRatio;
uniform vec4 standardColor;
uniform vec4 highlightColor;
in vec4 entryPoint; // = EntryPoint
out vec4 outColor;

void main(void)
{
    float stepSize = 1.0 / quality;
    vec3 exitPoint = ve_eyePosition.xyz;
    vec3 direction = exitPoint - entryPoint.xyz;
    vec3 deltaDir = normalize(direction) * stepSize;
    vec3 voxelCoord = entryPoint.xyz;
    vec4 accumulatedColor = vec4(0.0, 0.0, 0.0, 0.0);
    for(int i = 0; i < int(1.732 / stepSize); i++) { // 1.732 = cube diagonal
        if(any(lessThan(voxelCoord, vec3(0,0,0))) || any(greaterThan(voxelCoord, vec3(1,1,1)))) {
            break;
        }
        voxelCoord += deltaDir;
        float voxelValue = stepSize * texture(myTexture3D, voxelCoord).x;
        voxelValue = pow(voxelValue, contrast);
        voxelValue *= multiplier;
        float voxelValueClamped = clamp(voxelValue, 0.0, 1.0);
        float mixerValue = clamp(voxelValue*mixRatio, 0.0, 1.0);
        vec4 color = standardColor * mixerValue + highlightColor * (1 - mixerValue);
        accumulatedColor.a = color.a * voxelValueClamped + accumulatedColor.a * (1 - color.a * voxelValueClamped);
        accumulatedColor.rgb = (color.rgb * color.a * voxelValueClamped / accumulatedColor.a
                            + accumulatedColor.rgb * (1 - voxelValueClamped * color.a));
//        float newAlpha = color.a * voxelValueClamped;
//        accumulatedColor = newAlpha * color + (1 - newAlpha) * accumulatedColor;
        accumulatedColor = clamp(accumulatedColor, 0.0, 1.0);
    }
    outColor = accumulatedColor;
}
