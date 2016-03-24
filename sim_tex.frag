// Phong lighting in eye coordinates with texture and pure Phong specular.

// These are set by the .vert code, interpolated.
varying vec3 ec_vnormal, ec_vposition;

// This is set by the .c code.
uniform sampler2D mytexture;

void main()
{
    vec3 tcolor;
    float alpha;

    // perspective correction:
    tcolor = vec3(texture2D(mytexture,gl_TexCoord[0].st/gl_TexCoord[0].q));

    alpha = tcolor.x + tcolor.y + tcolor.z;
    if (alpha > 1.0f) { alpha = 1.0; }

    gl_FragColor = vec4(tcolor, alpha);
}