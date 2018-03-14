Shader "PBF/Thickness" {
	Properties{
		_ptlRadius("particle radius", Range(0.1, 1)) = 0.1
	}
	SubShader{
		Tags{ "Queue" = "Transparent" }
		Pass{

			ZWrite Off
			Blend One One

			Fog{ Mode off }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			//#pragma enable_d3d11_debug_symbols
			#pragma target 5.0

			#include "Lighting.cginc"
			#include "UnityCG.cginc"

			#pragma geometry geo

			struct Particle
			{
				float4 position;
				float4 velocity;
				float4 predictedPos;
				uint id;
				uint cellId;
				float pad;
				float pad2;
			};

			StructuredBuffer<Particle> _ParticlesBuffer;


			float _ptlRadius;

			struct v2g
			{
				float4	pos		: POSITION;
			};

			struct g2f
			{
				float4	pos		: POSITION;
				float2  tex0	: TEXCOORD0;
			};

			struct fout
			{
				float color		: COLOR;
			};

			v2g vert(uint vertex_id : SV_VertexID)
			{
				v2g o;

				float4 vPos = _ParticlesBuffer[vertex_id].position;
				o.pos = mul(unity_ObjectToWorld, vPos);

				return o;
			}

			[maxvertexcount(4)]
			void geo(point v2g p[1], inout TriangleStream<g2f> triStream)
			{
				float3 look = _WorldSpaceCameraPos - p[0].pos;
				look = normalize(look);

				//use camera's matrix
				float3 up = UNITY_MATRIX_IT_MV[1].xyz;

				float3 right = cross(up, look);

				float halfS = 0.5f * _ptlRadius;

				float4 v[4];
				v[0] = float4(p[0].pos + halfS * right - halfS * up, 1.0f);
				v[1] = float4(p[0].pos + halfS * right + halfS * up, 1.0f);
				v[2] = float4(p[0].pos - halfS * right - halfS * up, 1.0f);
				v[3] = float4(p[0].pos - halfS * right + halfS * up, 1.0f);


				float4x4 vp = UNITY_MATRIX_VP; //put the Quad into View Space
				g2f pIn;
				pIn.pos = mul(vp, v[0]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 0.0f));
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[1]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(1.0f, 1.0f));
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[2]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 0.0f));
				triStream.Append(pIn);

				pIn.pos = mul(vp, v[3]);
				pIn.tex0 = MultiplyUV(UNITY_MATRIX_TEXTURE0, float2(0.0f, 1.0f));
				triStream.Append(pIn);
			}

			fout frag(g2f input)
			{
				float gauConst = 0.3989422804014;

				fout OUT;

				float3 N;
				N.xy = input.tex0 * 2.0 - 1.0;
				float r2 = dot(N.xy, N.xy);
				if (r2 > 1.0) discard;	//kill pixels outside circle

				//Gaussian splat
				float dist = length(N.xy);
				float sigma = 4.0;
				float sigma2 = sigma * sigma;
				float gau = gauConst / sigma * exp(-dist*dist / 2.0*sigma2);

				/*
				if (gau < 0.1)
				{
					OUT.color = 0.5;
				}
				else
				{
					OUT.color = gau;
				}
				*/

				//OUT.color = float4(gau, gau, gau , 0.1);
				OUT.color = gau;

				return OUT;

			}
			ENDCG
		}
	}
}
