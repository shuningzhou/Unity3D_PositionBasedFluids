Shader "PBF/SphereBlurredDepth" {
	Properties{
		_MainTex("Base (RGB)", 2D) = "white" {}
		_DepthTex("Base (RGB)", 2D) = "white" {}
		scaleX("ScaleX", Float) = 0.1
		scaleY("ScaleY", Float) = 0.1
		//radius("Radius", Int) = 5
		//filterRadius("Filter Radius", Int) = 7
		//blurDepthFalloff("Falloff", Float) = 2.0
	}
	SubShader
	{
		Tags{ "RenderType" = "Opaque" }
		Pass
		{
			Cull Off
			ZWrite Off
			ZTest Always
			Fog{ Mode Off }

			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma fragment frag
			#pragma enable_d3d11_debug_symbols

			#include "UnityCG.cginc"
			sampler2D _DepthTex;
		    
			float blurDepthFalloff;
			float scaleX;
			float scaleY;
			int radius;
			int filterRadius;

			struct v2f
			{
				float4 pos		: POSITION;
				float2 coord	: TEXCOORD0;
			};

			struct fout
			{
				float color : COLOR;
			};

			v2f vert(appdata_base v)
			{
				v2f o;

				o.pos = UnityObjectToClipPos(v.vertex);
				o.coord = v.texcoord.xy;

				return o;
			}

			fout frag(v2f i)
			{
				fout OUT;
				float depth = tex2D(_DepthTex,  i.coord);

				float blurScale = 2.0 / radius;

				float sum = 0.0;
				float wsum = 0.0;	//sum of weights

				for (int x = -filterRadius; x <= filterRadius; x++)
				{
					float dSample = tex2D(_DepthTex,  i.coord + (float)x * float2(scaleX, scaleY));

					// spatial domain
					float r = (float)x * blurScale;
					float w = exp(-r*r);

					// range domain
					float r2 = (dSample - depth) * blurDepthFalloff;
					float g = exp(-r2*r2);

					sum += dSample * w * g;
					wsum += w * g;
				}

				if (wsum > 0.0)
				{
					sum /= wsum;
				}

				//	OUT.depth = sum;
				OUT.color = sum;
				//	sum = clamp(sum, 0.0, 0.999);
				//	OUT.color = EncodeFloatRGBA(sum);


				return OUT;
			}
			ENDCG
		}
	}
}


