package net.surina;

import net.surina.soundtouch.SoundTouch;
import net.surina.soundtouchexample.R;
import net.surina.soundtouchexample.R.id;
import net.surina.soundtouchexample.R.layout;
import android.app.Activity;
import android.os.Bundle;
import android.view.Menu;
import android.view.MenuItem;
import android.widget.TextView;

public class ExampleActivity extends Activity 
{
	TextView textViewConsole;
	StringBuilder consoleText = new StringBuilder();

	
	@Override
	protected void onCreate(Bundle savedInstanceState) 
	{
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_example);
		
		textViewConsole = (TextView)findViewById(R.id.textViewResult);
		
		// Check soundtouch
		checkLibVersion();
	}
	
	
	
	@Override
	protected void onDestroy()
	{
		textViewConsole = null;
	}
	
	
	/// Append text to console on the activity
	public void appendToConsole(String text)
	{
		if (textViewConsole == null) return;
		consoleText.append(text);
		consoleText.append("\n");
		textViewConsole.setText(consoleText);
	}
	
	
	
	/// print SoundTouch native library version onto console
	protected void checkLibVersion()
	{
		SoundTouch st = new SoundTouch();
		
		String ver = st.getVersionString();
		appendToConsole("SoundTouch native library version = " + ver);
	}

}
