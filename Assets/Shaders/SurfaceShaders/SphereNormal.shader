Shader "PBF/SphereNormal" {
	Properties{
		_DepthTex("Base (RGB)", 2D) = "white" {}
		textureWidth("Texture Width", Int) = 512
		textureHeight("Texture Height", Int) = 512
	}
	SubShader{
		Pass{
			Cull Off
			ZWrite Off
			ZTest Always
			Fog{ Mode Off }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma enable_d3d11_debug_symbols
			#pragma target 5.0

			#include "Lighting.cginc"
			#include "UnityCG.cginc"

			sampler2D _DepthTex;

			int textureWidth;
			int textureHeight;

			struct v2f
			{
				float4	pos			: POSITION;
				float2  coord		: TEXCOORD0;
			};

			struct fout
			{
				float4 color	: COLOR;
			};

			v2f vert(appdata_base v)
			{
				v2f o;

				o.pos = UnityObjectToClipPos(v.vertex);
				o.coord = v.texcoord.xy;

				return o;
			}

			float3 uvToEye(float2 texCoord, float z)
			{
				// Convert texture coordinate to homogeneous space
				float zFar = _ProjectionParams.z;
				float zNear = _ProjectionParams.y;

				float2 xyPos = (texCoord * 2.0 - 1.0);
				float a = zFar / (zFar - zNear);
				float b = zFar*zNear / (zNear - zFar);
				float rd = b / (z - a);
				return float3(xyPos.x, xyPos.y, -1.0) * rd;
			}

			fout frag(v2f input)
			{
				fout OUT;

				float2 uv = input.coord;
				float depth = tex2D(_DepthTex, input.coord);

				// reconstruct eye space pos from depth
				float3 eyePos = uvToEye(input.coord, depth);
				float ps = 0.5 / textureWidth;

				float2 tex = float2(uv.x + ps, uv.y);
				float2 tex2 = float2(uv.x - ps, uv.y);

				//calculate difference
				float3 ddx = uvToEye(tex, tex2D(_DepthTex, tex)) - eyePos;
				float3 ddx2 = eyePos - uvToEye(tex2, tex2D(_DepthTex, tex2));
				if (abs(ddx.z) > abs(ddx2.z))
				{
					ddx = ddx2;
				}

				ps = 0.5 / textureHeight;
				tex = float2(uv.x, uv.y + ps);
				tex2 = float2(uv.x, uv.y - ps);

				float3 ddy = uvToEye(tex, tex2D(_DepthTex, tex)) - eyePos;
				float3 ddy2 = eyePos - uvToEye(tex2, tex2D(_DepthTex, tex2));
				if (abs(ddy.z) > abs(ddy2.z))
				{
					ddy = ddy2;
				}

				//calculate normal
				float3 n = normalize(cross(ddx, ddy));

				//normal only
				OUT.color = float4(n, 1.0);

				return OUT;
			}
			ENDCG
		}
	}
}
